"""Shared mechanics for validated, all-or-nothing resource preparation."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from datetime import datetime, timezone
from typing import Any


class ResourcePreparationError(ValueError):
    """A scientific resource cannot be generated without invalid data or loss."""


def file_digest(path: Path, algorithm: str, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256(path: Path) -> str:
    return file_digest(path, "sha256")


def md5(path: Path) -> str:
    """Return an MD5 for compatibility with upstream scientific checksums."""
    return file_digest(path, "md5")


class JsonLinesRunLogger:
    """Append timestamped resource-preparation events to a canonical JSONL log."""

    def __init__(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        self.path = path

    def write(self, event: str, **details: Any) -> None:
        record = {
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "event": event,
            **details,
        }
        with self.path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(record, sort_keys=True) + "\n")

    def record_unfinished_previous_run(self) -> None:
        """Record recovery when the prior process ended before finalising its log."""
        if not self.path.is_file():
            return
        last_started: dict[str, Any] | None = None
        last_completed: dict[str, Any] | None = None
        with self.path.open("r", encoding="utf-8") as handle:
            for line_number, raw in enumerate(handle, 1):
                try:
                    record = json.loads(raw)
                except json.JSONDecodeError as exc:
                    raise ResourcePreparationError(
                        f"invalid JSON in canonical log {self.path} at line "
                        f"{line_number}"
                    ) from exc
                if record.get("event") == "run_started":
                    last_started = record
                    last_completed = None
                elif record.get("event") == "run_completed" and last_started is not None:
                    last_completed = record
        if last_started is not None and last_completed is None:
            self.write(
                "prior_run_interruption_detected",
                prior_run_started_utc=last_started.get("timestamp_utc"),
                recovery=(
                    "resuming validated complete files and byte-range partial downloads"
                ),
            )


def run_command(command: list[str], *, capture: bool = False) -> str:
    completed = subprocess.run(
        command,
        check=True,
        text=True,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.PIPE if capture else None,
    )
    return completed.stdout if capture else ""


def require_executable(name: str) -> str:
    path = shutil.which(name)
    if path is None:
        raise ResourcePreparationError(f"required executable is unavailable: {name}")
    return path


def resumable_download(
    *,
    curl: str,
    url: str,
    destination: Path,
    retries: int,
    partial_suffix: str,
    logger: JsonLinesRunLogger,
    expected_sha256: str | None = None,
) -> Path:
    """Download atomically and validate by configured SHA-256 or source size."""
    if expected_sha256 is not None and re.fullmatch(
        r"[0-9a-f]{64}", expected_sha256,
    ) is None:
        raise ResourcePreparationError("expected SHA-256 must be 64 lowercase hex digits")
    expected_bytes = None
    if expected_sha256 is None:
        headers = run_command(
            [curl, "--fail", "--location", "--silent", "--show-error", "--head", url],
            capture=True,
        )
        lengths = [
            int(line.split(":", 1)[1].strip())
            for line in headers.splitlines()
            if line.lower().startswith("content-length:")
            and line.split(":", 1)[1].strip().isdigit()
        ]
        if not lengths or lengths[-1] < 1:
            raise ResourcePreparationError(
                f"server did not provide a positive Content-Length: {url}"
            )
        expected_bytes = lengths[-1]
    if destination.is_file():
        if expected_sha256 is not None:
            observed_sha256 = sha256(destination)
            if observed_sha256 != expected_sha256:
                raise ResourcePreparationError(
                    f"cached download failed SHA-256 validation: {destination}"
                )
            logger.write(
                "download_reused", path=str(destination), sha256=observed_sha256,
            )
            return destination
        if destination.stat().st_size == expected_bytes:
            logger.write(
                "download_skipped",
                path=str(destination),
                reason="existing_size_matches_source",
                bytes=expected_bytes,
            )
            return destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(destination.name + partial_suffix)
    if destination.exists():
        if partial.exists():
            raise ResourcePreparationError(
                "both incomplete destination and partial file exist; inspect them "
                f"before retrying: {destination}, {partial}"
            )
        os.replace(destination, partial)
    if (
        partial.exists()
        and expected_bytes is not None
        and partial.stat().st_size > expected_bytes
    ):
        raise ResourcePreparationError(
            f"partial download is larger than the source ({partial.stat().st_size} > "
            f"{expected_bytes} bytes): {partial}"
        )
    if (
        expected_sha256 is not None
        and partial.is_file()
        and sha256(partial) == expected_sha256
    ):
        os.replace(partial, destination)
        logger.write(
            "download_completed",
            url=url,
            path=str(destination),
            sha256=expected_sha256,
        )
        return destination
    logger.write("download_started", url=url, path=str(destination))
    run_command([
        curl,
        "--fail",
        "--location",
        "--silent",
        "--show-error",
        "--retry",
        str(retries),
        "--retry-all-errors",
        "--continue-at",
        "-",
        "--output",
        str(partial),
        url,
    ])
    if expected_sha256 is not None:
        observed_sha256 = sha256(partial)
        if observed_sha256 != expected_sha256:
            raise ResourcePreparationError(
                f"download SHA-256 does not match for {url}: observed "
                f"{observed_sha256}, expected {expected_sha256}"
            )
    else:
        observed_bytes = partial.stat().st_size if partial.exists() else 0
        if observed_bytes != expected_bytes:
            raise ResourcePreparationError(
                f"download size does not match the source for {url}: observed "
                f"{observed_bytes}, expected {expected_bytes} bytes"
            )
    os.replace(partial, destination)
    details = (
        {"sha256": expected_sha256}
        if expected_sha256 is not None
        else {"bytes": expected_bytes}
    )
    logger.write("download_completed", url=url, path=str(destination), **details)
    return destination


def write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def required_file(value: str | Path, label: str) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_file() or path.stat().st_size == 0:
        raise ResourcePreparationError(f"{label} does not exist or is empty: {path}")
    return path


def create_staging_directory(output_directory: Path) -> Path:
    output = output_directory.expanduser().resolve()
    if output.exists():
        raise ResourcePreparationError(
            "output directory already exists; choose a new directory so an existing "
            f"scientific resource is never replaced implicitly: {output}"
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.", dir=output.parent))
    staging.chmod(output.parent.stat().st_mode & 0o777)
    return staging


def validate_output_filename(value: str, reserved_names: set[str]) -> str:
    if Path(value).name != value or not value.strip():
        raise ResourcePreparationError("output name must be a non-empty filename")
    if value in reserved_names:
        raise ResourcePreparationError(
            f"output name {value!r} is reserved for resource metadata"
        )
    return value


def unique_nonempty_values(values: list[str], label: str) -> set[str]:
    if not values or any(not value.strip() for value in values):
        raise ResourcePreparationError(f"{label} must contain non-empty values")
    if len(values) != len(set(values)):
        raise ResourcePreparationError(f"{label} must not contain duplicates")
    return set(values)


__all__ = [
    "ResourcePreparationError",
    "JsonLinesRunLogger",
    "create_staging_directory",
    "file_digest",
    "md5",
    "required_file",
    "require_executable",
    "resumable_download",
    "run_command",
    "sha256",
    "unique_nonempty_values",
    "validate_output_filename",
    "write_text",
]
