"""Safe path and executable resolution shared by PostGWAS modules."""

from __future__ import annotations

from pathlib import Path
import shutil
from typing import Type


def configured_output_path(
    root: str | Path,
    pattern: str,
    *,
    error_type: Type[Exception] = ValueError,
    **values,
) -> Path:
    """Render a relative configured pattern without allowing directory escape."""
    base = Path(root).expanduser().resolve()
    try:
        relative = Path(str(pattern).format(**values))
    except (KeyError, ValueError) as exc:
        raise error_type("Invalid output path pattern %r: %s" % (pattern, exc)) from exc
    if relative.is_absolute():
        raise error_type("Output path patterns must be relative: %s" % pattern)
    destination = (base / relative).resolve()
    if destination != base and base not in destination.parents:
        raise error_type("Output path pattern leaves its configured root: %s" % pattern)
    return destination


def validate_filename_component(
    value: object,
    label: str = "value",
    *,
    error_type: Type[Exception] = ValueError,
) -> str:
    """Validate a non-empty value used as one component of an output filename."""
    text = str(value).strip()
    if (
        not text
        or text in {".", ".."}
        or Path(text).name != text
        or "\x00" in text
    ):
        raise error_type(
            "%s must be a non-empty filename component without path separators"
            % label
        )
    return text


def expand_token_path(
    pattern: str,
    token: str,
    value: str,
    *,
    error_type: Type[Exception] = ValueError,
) -> Path:
    """Expand one required literal placeholder in a filesystem pattern."""
    if token not in pattern:
        raise error_type(
            "Path pattern must contain placeholder %r: %s" % (token, pattern)
        )
    return Path(pattern.replace(token, str(value))).expanduser().resolve()


def configured_output_matches(
    root: str | Path,
    pattern: str,
    *,
    error_type: Type[Exception] = ValueError,
    **values: object,
) -> list[Path]:
    """Find files using the same safe rendering as configured output paths."""
    rendered = configured_output_path(
        root, pattern, error_type=error_type, **values,
    )
    return list(rendered.parent.glob(rendered.name))


def resolve_executable(
    value: str | Path,
    label: str,
    *,
    error_type: Type[Exception] = RuntimeError,
) -> str:
    """Resolve an explicit executable path or search PATH for a command name."""
    candidate = Path(value).expanduser()
    if candidate.parent != Path(".") or candidate.is_absolute():
        resolved = candidate.resolve()
        if not resolved.is_file() or resolved.stat().st_size <= 0:
            raise error_type("%s does not exist or is empty: %s" % (label, resolved))
        return str(resolved)
    resolved = shutil.which(str(value))
    if not resolved:
        raise error_type("%s was not found: %s" % (label, value))
    return resolved


def require_nonempty_file(
    value: str | Path | None,
    label: str,
    *,
    error_type: Type[Exception] = ValueError,
) -> Path:
    """Resolve and require one regular, non-empty input file."""
    if value is None or not str(value).strip():
        raise error_type("%s is required" % label)
    path = Path(value).expanduser().resolve()
    if not path.is_file() or path.stat().st_size <= 0:
        raise error_type("%s does not exist or is empty: %s" % (label, path))
    return path


def remove_empty_directories(*values: str | Path) -> None:
    """Remove explicitly supplied directories only when they are empty."""
    for value in values:
        try:
            Path(value).rmdir()
        except (FileNotFoundError, OSError):
            pass


__all__ = [
    "configured_output_matches",
    "configured_output_path",
    "expand_token_path",
    "resolve_executable",
    "require_nonempty_file",
    "remove_empty_directories",
    "validate_filename_component",
]
