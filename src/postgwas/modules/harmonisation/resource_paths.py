"""Resolution and early validation of user-supplied external resource paths."""

from __future__ import annotations

import glob
from pathlib import Path
from string import Formatter

from postgwas.core.values import optional_text


_ALLOWED_TEMPLATE_FIELDS = frozenset(("build", "chromosome"))


def _template_glob(value: str) -> tuple[set[str], str]:
    """Return template fields and a safely escaped glob for their files."""
    fields: set[str] = set()
    parts: list[str] = []
    try:
        parsed = list(Formatter().parse(value))
    except ValueError as exc:
        raise ValueError("invalid path template %r: %s" % (value, exc)) from exc
    for literal, field, format_spec, conversion in parsed:
        parts.append(glob.escape(literal))
        if field is None:
            continue
        if field not in _ALLOWED_TEMPLATE_FIELDS:
            raise ValueError(
                "unsupported placeholder {%s}; allowed placeholders are "
                "{build} and {chromosome}" % field
            )
        if format_spec or conversion:
            raise ValueError(
                "placeholder {%s} must not use a format specifier or conversion"
                % field
            )
        fields.add(field)
        parts.append("*")
    return fields, "".join(parts)


def matching_external_resource_files(value: str | Path) -> tuple[Path, ...]:
    """Return exact or template-matched files without opening their contents."""
    text = str(value)
    fields, pattern = _template_glob(text)
    candidates = (
        (Path(match) for match in glob.glob(pattern))
        if fields
        else (Path(text),)
    )
    files = []
    for candidate in candidates:
        try:
            if candidate.is_file() and candidate.stat().st_size > 0:
                files.append(candidate.resolve())
        except OSError:
            continue
    return tuple(sorted(set(files), key=str))


def external_resource_template_fields(value: str | Path) -> frozenset[str]:
    """Return the validated placeholders used by one external-resource path."""
    fields, _pattern = _template_glob(str(value))
    return frozenset(fields)


def validate_external_resource_spec(value: str | Path) -> tuple[Path, ...]:
    """Require an existing file or an explicit template matching non-empty files."""
    text = str(value)
    fields, _pattern = _template_glob(text)
    files = matching_external_resource_files(text)
    if files:
        return files

    if fields:
        raise ValueError(
            "template matches no existing non-empty files: %s. Check the path "
            "and the {build}/{chromosome} placeholders." % text
        )

    prefix_matches = []
    for match in glob.glob(glob.escape(text) + "*"):
        candidate = Path(match)
        try:
            if candidate.is_file() and candidate.stat().st_size > 0:
                prefix_matches.append(candidate.resolve())
        except OSError:
            continue
    if prefix_matches:
        example = min(prefix_matches, key=lambda path: (len(str(path)), str(path)))
        raise ValueError(
            "the supplied value is a filename prefix, not a file: %s. Matching "
            "chromosome files exist, for example: %s. Enter the complete path "
            "template with the chromosome label replaced by {chromosome}, or "
            "enter one existing file that contains all chromosomes."
            % (text, example)
        )
    raise ValueError(
        "no existing non-empty file was found at %s. Enter one existing file, "
        "or an explicit per-chromosome template containing {chromosome}." % text
    )


def resolve_resource_file(
    input_file, grch_version, chromosome, *, must_exist=True,
):
    """Resolve one optional external-file template for a chromosome."""
    input_file = optional_text(input_file)
    if input_file is None:
        return None
    try:
        resolved = str(input_file).format(
            build=grch_version, chromosome=chromosome
        )
    except (KeyError, ValueError) as exc:
        raise ValueError(
            "Invalid external-file path template %r: %s" % (input_file, exc)
        ) from exc
    if must_exist and not Path(resolved).is_file():
        raise FileNotFoundError(
            "External file not found for chromosome %s: %s"
            % (chromosome, resolved)
        )
    return resolved


__all__ = [
    "external_resource_template_fields",
    "matching_external_resource_files",
    "resolve_resource_file",
    "validate_external_resource_spec",
]
