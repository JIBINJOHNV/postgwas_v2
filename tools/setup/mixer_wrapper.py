#!/usr/bin/env python3
"""Run an installed MiXeR entry point with its native library configured."""

from __future__ import annotations

import os
from pathlib import Path
import runpy
import sys
import tempfile


SUPPORTED_COMMANDS = {"mixer.py", "mixer_dev.py", "mixer_figures.py"}


def main() -> None:
    command = Path(sys.argv[0]).name
    if command not in SUPPORTED_COMMANDS:
        raise SystemExit("Unsupported MiXeR entry point: %s" % command)

    environment_prefix = Path(__file__).resolve().parents[1]
    mixer_environment = (
        environment_prefix / "share" / "postgwas" / "environments" / "mixer"
    )
    if sys.platform == "darwin" and mixer_environment.is_dir():
        mixer_python = (mixer_environment / "bin" / "python").resolve()
        current_python = Path(sys.executable).resolve()
        if not current_python.is_relative_to(mixer_environment.resolve()):
            if not mixer_python.is_file():
                raise SystemExit(
                    "MiXeR Intel-mac Python runtime is missing: %s"
                    % mixer_python
                )
            child_environment = os.environ.copy()
            child_environment["PYTHONNOUSERSITE"] = "1"
            child_environment.pop("PYTHONPATH", None)
            child_environment["DYLD_LIBRARY_PATH"] = str(
                mixer_environment / "lib"
            )
            child_environment.setdefault(
                "MPLCONFIGDIR",
                str(
                    Path(tempfile.gettempdir())
                    / ("postgwas-matplotlib-%s" % os.getuid())
                ),
            )
            os.execve(
                mixer_python,
                [str(mixer_python), str(Path(__file__).resolve()), *sys.argv[1:]],
                child_environment,
            )

    mixer_root = Path(
        os.environ.get(
            "MIXER_HOME",
            environment_prefix / "share" / "postgwas" / "mixer",
        )
    ).resolve()
    script = mixer_root / "precimed" / command
    library_name = "libbgmg.dylib" if sys.platform == "darwin" else "libbgmg.so"
    library = mixer_root / "lib" / library_name
    if not script.is_file():
        raise SystemExit("MiXeR script is missing: %s" % script)
    if not library.is_file():
        raise SystemExit("MiXeR native library is missing: %s" % library)

    os.environ["MIXER_HOME"] = str(mixer_root)
    os.environ["BGMG_SHARED_LIBRARY"] = str(library)
    sys.path.insert(0, str(script.parent))
    sys.argv[0] = str(script)
    runpy.run_path(str(script), run_name="__main__")


if __name__ == "__main__":
    main()
