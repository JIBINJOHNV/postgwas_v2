"""Fine-mapping logging that preserves detail without flooding the terminal."""

from __future__ import annotations

import logging
from contextlib import contextmanager
from pathlib import Path


@contextmanager
def detailed_file_logging(
    logger_name: str,
    log_file: str | Path,
    file_level: str,
    *,
    mode: str = "w",
):
    """Route one logger hierarchy to a detailed file and restore it afterward."""
    target = logging.getLogger(logger_name)
    previous_handlers = list(target.handlers)
    previous_level = target.level
    previous_propagate = target.propagate
    previous_disabled = target.disabled

    path = Path(log_file)
    path.parent.mkdir(parents=True, exist_ok=True)
    handler = logging.FileHandler(path, mode=mode, encoding="utf-8")
    handler.setLevel(str(file_level).upper())
    handler.setFormatter(
        logging.Formatter(
            "%(asctime)s | %(levelname)s | %(message)s",
        )
    )
    target.handlers = [handler]
    target.setLevel(str(file_level).upper())
    target.propagate = False
    target.disabled = False
    try:
        yield path
    finally:
        handler.flush()
        handler.close()
        target.handlers = previous_handlers
        target.setLevel(previous_level)
        target.propagate = previous_propagate
        target.disabled = previous_disabled


__all__ = ["detailed_file_logging"]
