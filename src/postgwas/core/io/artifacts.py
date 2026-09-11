"""Transactional publication for validated PostGWAS artifact sets."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence
from uuid import uuid4


def publish_artifact_set(
    pairs: Sequence[tuple[str | Path, str | Path]],
) -> None:
    """Publish validated temporary artifacts and restore prior outputs on error."""
    artifacts = [(Path(temporary), Path(final)) for temporary, final in pairs]
    if not artifacts:
        return
    temporary_paths = [temporary for temporary, _final in artifacts]
    final_paths = [final for _temporary, final in artifacts]
    if len(temporary_paths) != len(set(temporary_paths)):
        raise ValueError("Temporary artifact paths must be unique")
    if len(final_paths) != len(set(final_paths)):
        raise ValueError("Final artifact paths must be unique")
    missing = [str(path) for path in temporary_paths if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "Validated temporary artifacts are missing: %s" % ", ".join(missing)
        )

    token = uuid4().hex
    backups = {
        final: final.with_name(".%s.%s.backup" % (final.name, token))
        for final in final_paths
        if final.exists()
    }
    published: list[Path] = []
    try:
        for final, backup in backups.items():
            final.replace(backup)
        for temporary, final in artifacts:
            temporary.replace(final)
            published.append(final)
    except OSError as exc:
        rollback_errors = []
        for final in reversed(published):
            try:
                final.unlink(missing_ok=True)
            except OSError as rollback_exc:
                rollback_errors.append(str(rollback_exc))
        for final, backup in backups.items():
            if not backup.exists():
                continue
            try:
                backup.replace(final)
            except OSError as rollback_exc:
                rollback_errors.append(str(rollback_exc))
        suffix = (
            " Rollback also failed: %s" % "; ".join(rollback_errors)
            if rollback_errors
            else " Existing outputs were restored."
        )
        raise RuntimeError(
            "Cannot publish the validated artifact set: %s.%s" % (exc, suffix)
        ) from exc
    for backup in backups.values():
        backup.unlink(missing_ok=True)


__all__ = ["publish_artifact_set"]
