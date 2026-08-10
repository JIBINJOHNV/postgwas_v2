#!/usr/bin/env python3
"""Download and verify the FinnGen R13 annotated-variant resource."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import sys
from datetime import datetime, timezone


# Edit resource and download settings here.
SOURCE_URL = (
    "https://storage.googleapis.com/finngen-public-data-r13/annotations/"
    "finngen_R13_annotated_variants_v0.gz"
)
OUTPUT_DIRECTORY = Path(
    "/Users/JJOHN41/Documents/software_resources/resourses/postgwas_v2/"
    "finngen_r13_annotations"
)
OUTPUT_FILENAME = "finngen_R13_annotated_variants_v0.gz"
PARTIAL_SUFFIX = ".partial"
WGET_EXECUTABLE = "/Users/JJOHN41/miniconda3/bin/wget"
DOWNLOAD_RETRIES = 20
NETWORK_TIMEOUT_SECONDS = 60

# Values published in the Google Cloud Storage response headers.
EXPECTED_BYTES = 29_768_495_399
EXPECTED_MD5_HEX = "367d21b33fdda7adb52786313ced0eab"


def timestamp() -> str:
    """Return an ISO-8601 UTC timestamp for concise progress messages."""
    return datetime.now(timezone.utc).isoformat()


def md5_hex(path: Path) -> str:
    """Calculate a file MD5 without loading the large resource into memory."""
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_download(path: Path) -> None:
    """Validate the exact published object size and MD5 checksum."""
    actual_bytes = path.stat().st_size
    if actual_bytes != EXPECTED_BYTES:
        raise RuntimeError(
            f"incorrect file size for {path}: expected {EXPECTED_BYTES}, "
            f"found {actual_bytes}"
        )
    print(f"[{timestamp()}] Size validated: {actual_bytes:,} bytes", flush=True)

    actual_md5 = md5_hex(path)
    if actual_md5 != EXPECTED_MD5_HEX:
        raise RuntimeError(
            f"MD5 mismatch for {path}: expected {EXPECTED_MD5_HEX}, "
            f"found {actual_md5}"
        )
    print(f"[{timestamp()}] MD5 validated: {actual_md5}", flush=True)


def main() -> int:
    """Resume the download, validate it, and publish it atomically."""
    wget = shutil.which(WGET_EXECUTABLE)
    if wget is None:
        raise RuntimeError(f"wget executable was not found: {WGET_EXECUTABLE}")

    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    destination = OUTPUT_DIRECTORY / OUTPUT_FILENAME
    partial = Path(f"{destination}{PARTIAL_SUFFIX}")

    if destination.exists():
        print(f"[{timestamp()}] Existing file found: {destination}", flush=True)
        validate_download(destination)
        print(f"[{timestamp()}] Download already complete", flush=True)
        return 0

    if partial.exists() and partial.stat().st_size > EXPECTED_BYTES:
        raise RuntimeError(
            f"partial file is larger than the official object: {partial}"
        )

    print(
        f"[{timestamp()}] Downloading {SOURCE_URL}\n"
        f"[{timestamp()}] Destination: {destination}\n"
        f"[{timestamp()}] Existing partial bytes: "
        f"{partial.stat().st_size if partial.exists() else 0:,}",
        flush=True,
    )
    subprocess.run(
        [
            wget,
            "--no-verbose",
            "--continue",
            "--tries",
            str(DOWNLOAD_RETRIES),
            "--timeout",
            str(NETWORK_TIMEOUT_SECONDS),
            "--output-document",
            str(partial),
            SOURCE_URL,
        ],
        check=True,
    )

    validate_download(partial)
    os.replace(partial, destination)
    print(f"[{timestamp()}] Download complete: {destination}", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"[{timestamp()}] ERROR: {error}", file=sys.stderr, flush=True)
        raise SystemExit(1) from error
