"""Content-verified reuse of installer components, without importing PostGWAS.

Bootstrap runs before the package is available. Receipts describe software,
not scientific results, and are published only after the installer verifies
the complete selected software stack.
"""

from __future__ import annotations

import argparse
import glob
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import sys
import tempfile


def snapshot(patterns: list[str]) -> dict:
    result = {}
    for pattern in patterns:
        roots = sorted(glob.glob(pattern))
        if not roots:
            raise OSError(f"Required component path is missing: {pattern}")
        for root_text in roots:
            root = Path(root_text).absolute()
            paths = [root, *sorted(root.rglob("*"))] if root.is_dir() else [root]
            for path in paths:
                # Interpreter caches do not define installed source/software.
                if "__pycache__" in path.parts or path.suffix in {".pyc", ".pyo"}:
                    continue
                entry = {"mode": path.stat().st_mode & 0o777}
                if path.is_symlink():
                    entry["link"] = os.readlink(path)
                if path.is_file():
                    digest = hashlib.sha256()
                    with path.open("rb") as handle:
                        for block in iter(lambda: handle.read(1024 * 1024), b""):
                            digest.update(block)
                    entry["sha256"] = digest.hexdigest()
                elif not path.is_dir():
                    raise OSError(f"Unsupported component path: {path}")
                result[str(path)] = entry
    return result


def atomic_write(path: Path, document: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(dir=path.parent, prefix=".receipt-")
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(document, handle, sort_keys=True, indent=2)
            handle.write("\n")
        os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


def component_document(args) -> dict:
    if not args.input or not args.output:
        raise ValueError("A component requires declared inputs and outputs")
    return {
        "format": 1,
        "prefix": str(args.prefix.resolve()),
        "platform": [platform.system(), platform.machine()],
        "values": args.value,
        "input_paths": args.input,
        "output_paths": args.output,
        "inputs": snapshot(args.input),
        "outputs": snapshot(args.output),
    }


def commit_receipts(prefix: Path, staging: Path) -> None:
    pending = []
    for path in sorted(staging.glob("*.json")):
        document = json.loads(path.read_text(encoding="utf-8"))
        if (document["prefix"] != str(prefix.resolve())
                or document["inputs"] != snapshot(document["input_paths"])
                or document["outputs"] != snapshot(document["output_paths"])):
            raise ValueError(f"Component changed during installation: {path.stem}")
        pending.append((path.name, document))
    for name, document in pending:
        atomic_write(prefix / "share/postgwas/install-state" / name, document)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("operation", choices=("check", "record", "commit"))
    parser.add_argument("--prefix", type=Path, default=Path(sys.prefix))
    parser.add_argument("--staging", type=Path, required=True)
    parser.add_argument("--component")
    parser.add_argument("--input", action="append", default=[])
    parser.add_argument("--output", action="append", default=[])
    parser.add_argument("--value", action="append", default=[])
    args = parser.parse_args()
    if args.operation != "commit" and not re.fullmatch(r"[a-z][a-z0-9-]*", args.component or ""):
        parser.error("--component must be a lowercase component name")
    try:
        if args.operation == "commit":
            commit_receipts(args.prefix, args.staging)
            return 0
        document = component_document(args)
        name = args.component + ".json"
        if args.operation == "record":
            atomic_write(args.staging / name, document)
        else:
            installed = args.prefix / "share/postgwas/install-state" / name
            if json.loads(installed.read_text(encoding="utf-8")) != document:
                print(f"Updating {args.component}: recipe, dependencies or installed files changed.")
                return 1
            print(f"Reusing verified component: {args.component}")
    except (OSError, ValueError, KeyError, TypeError) as exc:
        if args.operation == "check":
            print(f"Updating {args.component}: no matching complete build record ({exc}).")
            return 1
        print(f"ERROR: installation state {args.operation} failed: {exc}")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
