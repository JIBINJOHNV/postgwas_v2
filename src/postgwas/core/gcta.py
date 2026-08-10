"""Shared GCTA executable validation used by all GCTA-backed modules."""

from __future__ import annotations

import re
import subprocess
import tempfile
from pathlib import Path
from typing import Type


def _version_tuple(value: str) -> tuple[int, ...]:
    return tuple(int(part) for part in value.split("."))


def require_supported_gcta(
    executable: str,
    module_config,
    logger,
    timeout_seconds: float,
    *,
    error_type: Type[Exception] = RuntimeError,
    configuration_path: str,
    analysis_name: str,
) -> str:
    """Read the GCTA banner and enforce a module's configured minimum version.

    GCTA has no dedicated version-only mode: supported releases may print a
    valid banner and still return nonzero for the configured probe argument.
    Probe output is isolated in a temporary directory, and the configured
    banner pattern is authoritative. Analysis commands still require zero.
    """
    try:
        with tempfile.TemporaryDirectory(
            prefix=module_config.version_probe_temporary_prefix,
        ) as temporary_directory:
            output_prefix = (
                Path(temporary_directory) / module_config.version_probe_output_name
            )
            command = [
                executable,
                *module_config.version_arguments,
                module_config.version_probe_output_argument,
                str(output_prefix),
            ]
            logger.record(
                "INPUT",
                "external_command",
                purpose="Read GCTA version",
                command=command,
            )
            result = subprocess.run(
                command,
                capture_output=True,
                text=True,
                timeout=timeout_seconds,
            )
    except subprocess.TimeoutExpired as exc:
        raise error_type(
            "Read GCTA version exceeded the configured timeout of %s seconds"
            % timeout_seconds
        ) from exc
    except OSError as exc:
        raise error_type(
            "Read GCTA version could not create temporary output or start: %s" % exc
        ) from exc
    output = "\n".join(value for value in (result.stdout, result.stderr) if value)
    match = re.search(module_config.version_pattern, output)
    if match is None:
        raise error_type(
            "GCTA version probe exited with status %d and its output did not match "
            "%s.version_pattern: %r"
            % (result.returncode, configuration_path, output.strip())
        )
    version = match.group(1)
    observed = _version_tuple(version)
    minimum = _version_tuple(module_config.minimum_gcta_version)
    width = max(len(observed), len(minimum))
    if observed + (0,) * (width - len(observed)) < minimum + (0,) * (
        width - len(minimum)
    ):
        raise error_type(
            "GCTA %s is too old; %s requires at least %s."
            % (version, analysis_name, module_config.minimum_gcta_version)
        )
    logger.record(
        "OBSERVED", "gcta_version", version=version, executable=executable,
        probe_exit_code=result.returncode,
    )
    return version


__all__ = ["require_supported_gcta"]
