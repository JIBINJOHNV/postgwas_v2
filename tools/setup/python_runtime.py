"""Install canonical pip requirements and exercise installed public commands.

Called only through the selected environment's isolated Python by the installer.
The environment YAML remains the sole source of pip package/version choices.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys

import yaml


def install_requirements(environment_file: Path) -> None:
    configuration = yaml.safe_load(environment_file.read_text(encoding="utf-8"))
    dependencies = configuration.get("dependencies") if isinstance(configuration, dict) else None
    if not isinstance(dependencies, list):
        raise ValueError("Environment YAML must contain a dependencies list")
    groups = [entry["pip"] for entry in dependencies if isinstance(entry, dict) and "pip" in entry]
    if len(groups) != 1 or not isinstance(groups[0], list) or not groups[0]:
        raise ValueError("Environment YAML must contain one nonempty pip requirements list")
    requirements = groups[0]
    if any(not isinstance(value, str) or not value.strip() or value.startswith("-")
           for value in requirements):
        raise ValueError("Each pip requirement must be a nonempty package requirement")
    subprocess.run(
        [sys.executable, "-I", "-m", "pip", "--isolated", "install", *requirements],
        check=True,
    )


def verify_commands() -> None:
    from postgwas.pipeline.registry import REGISTRY

    failures = []
    # --validate is a public dispatcher branch outside the module registry.
    commands = [*REGISTRY.commands(), "--validate"]
    for command in commands:
        result = subprocess.run(
            [sys.executable, "-I", "-m", "postgwas", command, "--help"],
            capture_output=True, text=True, check=False,
        )
        if result.returncode:
            failures.append(command)
            print(f"ERROR: postgwas {command} --help failed (exit {result.returncode})",
                  file=sys.stderr)
            print(result.stdout + result.stderr, file=sys.stderr)
    if failures:
        raise RuntimeError("Python command verification failed: " + ", ".join(failures))
    print(f"Verified {len(commands)} installed public command entry points.")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("operation", choices=("install", "verify"))
    parser.add_argument("--environment-file", type=Path)
    args = parser.parse_args()
    if args.operation == "install" and args.environment_file is None:
        parser.error("install requires --environment-file")
    try:
        if args.operation == "install":
            install_requirements(args.environment_file)
        else:
            verify_commands()
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as exc:
        print(f"ERROR: Python runtime {args.operation} failed: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
