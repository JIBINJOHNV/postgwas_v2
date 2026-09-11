"""Small read-only contracts for configured reference-resource files.

Availability checks do not establish scientific validity. File identities and
directory membership are checked again on use; parsed static content is reused
only inside the active input-validation session.
"""

from __future__ import annotations

import csv
from pathlib import Path

import yaml

from postgwas.core.input_validation import record_file_validation, validate_once
from postgwas.core.io.delimiters import open_text
from postgwas.core.paths import require_nonempty_file


def require_file_inventory(paths, label, *, missing_message, error_type=ValueError):
    """Require all configured files, reporting every independently absent file."""
    resolved = []
    missing = []
    for value in paths:
        try:
            resolved.append(require_nonempty_file(value, label, error_type=error_type))
        except error_type:
            missing.append(str(value))
    if missing:
        raise error_type("%s: %s" % (missing_message, ", ".join(missing)))
    return tuple(resolved)


def contained_resource_path(root, relative, *, label, root_label, error_type=ValueError):
    """Resolve a configured resource without allowing it to escape its root."""
    root = Path(root).expanduser().resolve()
    path = (root / relative).resolve()
    if path != root and root not in path.parents:
        message = "%s leaves %s" % (label, root_label)
        record_file_validation(
            path, label, checks=("configured root containment",),
            status="failed", message=message,
        )
        raise error_type(message)
    return path


def validate_resource_directories(root, directories, *, label, bundle, error_type=ValueError):
    """Inventory nonempty files in each required directory, without parsing them.

    Re-enumeration is intentional: a list of previously known file identities
    cannot detect a newly added resource or a newly empty directory.
    """
    root = Path(root).expanduser().resolve()
    missing, empty = [], []
    validated_directories, inventory = {}, {}
    for relative in directories:
        directory = root / relative
        if not directory.is_dir():
            missing.append(relative)
            record_file_validation(
                directory, label, checks=("directory availability",),
                status="failed", message="Configured directory does not exist.",
            )
            continue
        files = sorted(entry.resolve() for entry in directory.rglob("*") if entry.is_file())
        if not files:
            empty.append(relative)
            record_file_validation(
                directory, label, checks=("nonempty resource inventory",),
                status="failed", message="Configured directory contains no files.",
            )
            continue
        validated_directories[relative] = directory.resolve()
        for path in files:
            contained_resource_path(
                root, path, label=label, root_label="the configured resource directory",
                error_type=error_type,
            )
            validated = require_nonempty_file(path, label, error_type=error_type)
            inventory[str(validated.relative_to(root))] = validated
    if missing:
        raise error_type(
            "%s is missing configured directories: %s" % (bundle, ", ".join(missing))
        )
    if empty:
        raise error_type(
            "%s has configured resource directories without any files: %s"
            % (bundle, ", ".join(empty))
        )
    return validated_directories, inventory


def read_yaml_manifest(path, model_type, label, *, error_type=ValueError):
    """Parse a YAML mapping and its explicitly supplied schema once per run."""
    path = require_nonempty_file(path, label, error_type=error_type)
    checks = ("YAML mapping", "configured manifest schema")

    def inspect():
        try:
            try:
                document = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
            except (OSError, UnicodeError, yaml.YAMLError) as exc:
                raise error_type("Cannot read %s %s: %s" % (label, path, exc)) from exc
            if not isinstance(document, dict):
                raise error_type("%s must be a YAML mapping: %s" % (label, path))
            try:
                result = model_type.model_validate(document)
            except Exception as exc:
                raise error_type("Invalid %s %s: %s" % (label, path, exc)) from exc
        except error_type as exc:
            record_file_validation(path, label, checks=checks, status="failed", message=str(exc))
            raise
        record_file_validation(
            path, label, checks=checks, metrics={"declared_fields": len(document)},
            message="Declarations validated; referenced file contents are not established by the manifest.",
        )
        return result

    result = validate_once((path,), {
        "validator": "yaml_manifest",
        "model": "%s.%s" % (model_type.__module__, model_type.__qualname__),
        "schema": model_type.model_json_schema(),
    }, inspect, error_type=error_type)
    return result.model_copy(deep=True)


def validate_table_header(path, *, delimiter, required_columns, label, error_type=ValueError):
    """Check named columns and the presence of one record, not all row values."""
    path = require_nonempty_file(path, label, error_type=error_type)
    required_columns = tuple(required_columns)
    checks = ("configured header columns", "at least one data record")

    def inspect():
        try:
            with open_text(path) as handle:
                reader = csv.DictReader(handle, delimiter=delimiter)
                header = reader.fieldnames or []
                missing = [name for name in required_columns if name not in header]
                if missing:
                    raise error_type(
                        "%s is missing configured columns: %s" % (label, ", ".join(missing))
                    )
                if next(reader, None) is None:
                    raise error_type("%s contains no records: %s" % (label, path))
        except (OSError, UnicodeError, csv.Error) as exc:
            message = "Cannot read %s %s: %s" % (label, path, exc)
            record_file_validation(path, label, checks=checks, status="failed", message=message)
            raise error_type(message) from exc
        except error_type as exc:
            record_file_validation(path, label, checks=checks, status="failed", message=str(exc))
            raise
        record_file_validation(
            path, label, checks=checks, metrics={"columns": len(header)},
            message="Header and first-record presence only; remaining records are not validated.",
        )
        return tuple(header)

    return validate_once((path,), {
        "validator": "delimited_header_and_first_record",
        "delimiter": delimiter, "required_columns": required_columns,
    }, inspect, error_type=error_type)


def read_unique_line_names(path, label, *, error_type=ValueError):
    """Read unique stripped nonblank lines; spaces and '#' remain name content."""
    path = require_nonempty_file(path, label, error_type=error_type)
    checks = ("nonempty unique line names",)

    def inspect():
        try:
            names = tuple(
                line.strip() for line in path.read_text(encoding="utf-8").splitlines()
                if line.strip()
            )
            if not names or len(names) != len(set(names)):
                raise error_type("%s must contain unique non-empty names: %s" % (label, path))
        except (OSError, UnicodeError) as exc:
            message = "Cannot read %s %s: %s" % (label, path, exc)
            record_file_validation(path, label, checks=checks, status="failed", message=message)
            raise error_type(message) from exc
        except error_type as exc:
            record_file_validation(path, label, checks=checks, status="failed", message=str(exc))
            raise
        record_file_validation(path, label, checks=checks, metrics={"names": len(names)})
        return names

    return validate_once(
        (path,), {"validator": "unique_nonblank_line_names"}, inspect,
        error_type=error_type,
    )
