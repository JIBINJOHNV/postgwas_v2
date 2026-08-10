#!/usr/bin/env python3
"""Download and verify the latest NCBI dbSNP ALFA frequency VCF."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone


# Edit download and resource settings here.
SOURCE_DIRECTORY_URL = (
    "https://ftp.ncbi.nih.gov/snp/population_frequency/latest_release"
)
OUTPUT_DIRECTORY = Path(
    "/Users/JJOHN41/Documents/software_resources/resourses/postgwas_v2/"
    "dbsnp_alfa_latest"
)
DATA_FILENAME = "freq.vcf.gz"
INDEX_FILENAME = "freq.vcf.gz.tbi"
CHECKSUM_SUFFIX = ".md5"
PARTIAL_SUFFIX = ".partial"
WGET_EXECUTABLE = "/Users/JJOHN41/miniconda3/bin/wget"
DOWNLOAD_RETRIES = 20
NETWORK_TIMEOUT_SECONDS = 60
HASH_CHUNK_BYTES = 8 * 1024 * 1024

# Current release metadata, used to detect a changed `latest_release` during a run.
EXPECTED_DATA_BYTES = 19_721_864_743
EXPECTED_INDEX_BYTES = 2_971_690


def timestamp() -> str:
    """Return an ISO-8601 UTC timestamp."""
    return datetime.now(timezone.utc).isoformat()


def md5_hex(path: Path) -> str:
    """Calculate MD5 without loading a large file into memory."""
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(HASH_CHUNK_BYTES), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(wget: str, filename: str, destination: Path) -> None:
    """Resume one download into a partial file and publish it atomically."""
    partial = Path(f"{destination}{PARTIAL_SUFFIX}")
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
            f"{SOURCE_DIRECTORY_URL}/{filename}",
        ],
        check=True,
    )
    os.replace(partial, destination)


def read_published_md5(path: Path, expected_filename: str) -> str:
    """Parse an NCBI MD5 sidecar and verify that it names the expected file."""
    fields = path.read_text(encoding="utf-8").strip().split()
    if len(fields) != 2 or fields[1].lstrip("*") != expected_filename:
        raise RuntimeError(f"invalid MD5 sidecar content: {path}")
    if re.fullmatch(r"[0-9a-fA-F]{32}", fields[0]) is None:
        raise RuntimeError(f"invalid MD5 digest in {path}")
    return fields[0].lower()


def validate(path: Path, expected_bytes: int, expected_md5: str) -> None:
    """Validate exact size and NCBI-published checksum."""
    actual_bytes = path.stat().st_size
    if actual_bytes != expected_bytes:
        raise RuntimeError(
            f"incorrect size for {path}: expected {expected_bytes:,}, "
            f"found {actual_bytes:,}"
        )
    actual_md5 = md5_hex(path)
    if actual_md5 != expected_md5:
        raise RuntimeError(
            f"MD5 mismatch for {path}: expected {expected_md5}, found {actual_md5}"
        )
    print(
        f"[{timestamp()}] Validated {path.name}: {actual_bytes:,} bytes, {actual_md5}",
        flush=True,
    )


def main() -> int:
    """Download the ALFA VCF, index, and checksums, then validate both resources."""
    wget = shutil.which(WGET_EXECUTABLE)
    if wget is None:
        raise RuntimeError(f"wget executable was not found: {WGET_EXECUTABLE}")

    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    resources = (
        (DATA_FILENAME, EXPECTED_DATA_BYTES),
        (INDEX_FILENAME, EXPECTED_INDEX_BYTES),
    )

    for filename, _ in resources:
        checksum_filename = f"{filename}{CHECKSUM_SUFFIX}"
        checksum_path = OUTPUT_DIRECTORY / checksum_filename
        print(f"[{timestamp()}] Fetching {checksum_filename}", flush=True)
        download(wget, checksum_filename, checksum_path)

    for filename, expected_bytes in resources:
        destination = OUTPUT_DIRECTORY / filename
        checksum_path = OUTPUT_DIRECTORY / f"{filename}{CHECKSUM_SUFFIX}"
        expected_md5 = read_published_md5(checksum_path, filename)

        if destination.exists():
            print(f"[{timestamp()}] Checking existing {destination}", flush=True)
        else:
            partial = Path(f"{destination}{PARTIAL_SUFFIX}")
            partial_bytes = partial.stat().st_size if partial.exists() else 0
            if partial_bytes > expected_bytes:
                raise RuntimeError(
                    f"partial file is larger than the current release: {partial}"
                )
            print(
                f"[{timestamp()}] Downloading {filename}; "
                f"existing partial bytes: {partial_bytes:,}",
                flush=True,
            )
            download(wget, filename, destination)

        validate(destination, expected_bytes, expected_md5)

    print(
        f"[{timestamp()}] NCBI dbSNP ALFA download complete: {OUTPUT_DIRECTORY}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"[{timestamp()}] ERROR: {error}", file=sys.stderr, flush=True)
        raise SystemExit(1) from error
