"""Discover and validate an R installation without mutating the parent environment."""

import json
import os
from pathlib import Path
import subprocess

from postgwas.core.paths import resolve_executable


def resolve_r_runtime(rscript, timeout_seconds, required_packages, *, label):
    """Prefer the executable's own library and validate its required namespaces.

    R's .Library is the installation's native library. Its precedence prevents
    incompatible user libraries from shadowing bundled packages; discovered
    libraries remain available afterwards. See R's base::.libPaths manual.
    The returned complete environment is for the child process only.
    """
    resolved_rscript = resolve_executable(rscript, "Rscript executable")
    discovery = subprocess.run(
        [resolved_rscript, "--vanilla", "-e",
         "cat(normalizePath(.Library, mustWork=TRUE), '\\n', sep=''); "
         "cat(paste(.libPaths(), collapse=.Platform$path.sep), '\\n', sep='')"],
        check=False, capture_output=True, text=True, timeout=timeout_seconds,
    )
    if discovery.returncode != 0:
        detail = (discovery.stderr or discovery.stdout).strip()
        raise RuntimeError(
            "Configured Rscript could not report its library paths: "
            f"{resolved_rscript}; {detail}"
        )
    lines = [line.strip() for line in discovery.stdout.splitlines() if line.strip()]
    if len(lines) < 2:
        raise RuntimeError(
            "Configured Rscript returned incomplete library-path metadata: "
            f"{resolved_rscript}"
        )
    ordered_libraries = list(dict.fromkeys([lines[0], *lines[1].split(os.pathsep)]))
    environment = os.environ.copy()
    # R_LIBS precedes R_LIBS_USER at startup; set both for consistent ordering.
    environment["R_LIBS"] = os.pathsep.join(ordered_libraries)
    environment["R_LIBS_USER"] = os.pathsep.join(ordered_libraries)
    environment["PATH"] = os.pathsep.join(
        [str(Path(resolved_rscript).parent), environment.get("PATH", "")]
    ).rstrip(os.pathsep)
    package_vector = ", ".join(json.dumps(package) for package in required_packages)
    validation = subprocess.run(
        [resolved_rscript, "--vanilla", "-e",
         f"required <- c({package_vector}); "
         "failed <- required[!vapply(required, function(package) "
         "isTRUE(requireNamespace(package, quietly=TRUE)), logical(1))]; "
         "if (length(failed)) stop(paste('unavailable packages:', "
         "paste(failed, collapse=', '))); "
         "cat(R.version.string, '\\n', sep=''); "
         "for (package in required) cat(package, as.character(packageVersion(package)), '\\n', sep='\\t')"],
        check=False, capture_output=True, text=True, timeout=timeout_seconds,
        env=environment,
    )
    if validation.returncode != 0:
        detail = (validation.stderr or validation.stdout).strip()
        raise RuntimeError(
            f"Configured R/{label} runtime failed dependency validation: "
            f"{resolved_rscript}; {detail}"
        )
    metadata = validation.stdout.strip().splitlines()
    return {
        "rscript": resolved_rscript,
        "version": metadata[0] if metadata else "",
        "package_versions": dict(line.strip().split("\t", 1) for line in metadata[1:]),
        "library_paths": ordered_libraries,
        "environment": environment,
    }
