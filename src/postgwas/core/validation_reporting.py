"""File and explicit bundle presentation of checks, without reopening inputs.

The audit combines distinct check observations under each file. Grouping never upgrades their
scope: an availability check is not content validation, and a failed companion
or cross-file requirement cannot be hidden by a successful file-existence check.
"""

from __future__ import annotations

from collections import Counter
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, field, replace
import json
import os
from pathlib import Path
import re

from postgwas.core.input_validation import (
    InputValidationSession, current_validation_session, current_validation_recorder,
    record_file_validation,
)
from postgwas.core.io.reports import write_yaml_report
from postgwas.core.paths import configured_output_path, validate_filename_component
from postgwas.core.ui import print_screen_block
from postgwas.core.ui.screen import SUMMARY_CARD_LABEL_WIDTH, screen_field, screen_line


# Serialization/check vocabulary and standard companion extensions, not
# scientific policies or filename-based guesses about the parent file's format.
_STATUS_ORDER = {"failed": 0, "blocked": 1, "warning": 2, "deferred": 3, "passed": 4}
_STATUS_KIND = {"failed": "error", "blocked": "error", "warning": "warning",
                "deferred": "decision", "passed": "success"}
_AVAILABILITY = frozenset({
    "regular file", "non-empty file", "regular nonempty file", "readability",
    "availability", "directory availability", "nonempty resource inventory",
    "directory exists and contains a file", "unchanged since validation",
})
_INDEX_SUFFIXES = (".tbi", ".csi", ".fai", ".gzi")
_DIRECT_REPORT_KIND = "postgwas.direct.input_validation"
_DIRECT_DISPLAY: ContextVar[object | None] = ContextVar(
    "postgwas_direct_file_validation_display", default=None,
)


@dataclass
class FileValidationGroup:
    path: str
    records: list = field(default_factory=list)
    companions: dict[str, list] = field(default_factory=dict)

    @property
    def all_records(self):
        return [*self.records, *(record for items in self.companions.values() for record in items)]

    @property
    def status(self):
        return min((record.status for record in self.all_records), key=_STATUS_ORDER.__getitem__)


@dataclass(frozen=True)
class FileValidationBundle:
    """Presentation-only grouping for an explicitly declared file inventory."""

    role: str
    paths: tuple[str, ...]
    fields: tuple[tuple[str, str, object, bool], ...]
    covered_checks: frozenset[str]
    covered_metric_keys: frozenset[str]
    availability_only: bool


def group_file_validations(records):
    """Group exact recorded paths; do not stat, read, normalise alleles, or infer format."""
    by_path = {}
    requirements = []
    for record in records:
        paths = (record.path,) if record.path is not None else record.metrics.get("paths", ())
        if not paths:
            requirements.append(record)
        for path in dict.fromkeys(paths):
            by_path.setdefault(str(path), []).append(record)
    groups = {}
    for path, items in by_path.items():
        parent = next((path[:-len(suffix)] for suffix in _INDEX_SUFFIXES
                       if path.endswith(suffix) and path[:-len(suffix)] in by_path), path)
        group = groups.setdefault(parent, FileValidationGroup(parent))
        if parent == path:
            group.records.extend(items)
        else:
            group.companions[path] = items
    # Stable input order on success; problems first. No warning/failure is lost
    # if the user limits the number of successful file sections on screen.
    return sorted(groups.values(), key=lambda group: _STATUS_ORDER[group.status]), requirements


def _value_text(value, limit):
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, int):
        return f"{value:,}"
    if isinstance(value, dict):
        items = ["%s=%s" % (key, _value_text(item, limit)) for key, item in value.items()]
    elif isinstance(value, (list, tuple)):
        items = [str(item) for item in value]
    else:
        return str(value)
    shown = ", ".join(items[:limit]) or "none"
    return shown if len(items) <= limit else "%s (+%d more; see audit)" % (shown, len(items) - limit)


def _contig_text(values, limit):
    """Compress consecutive labels for display without chromosome aliases or inference."""
    if not isinstance(values, (list, tuple)):
        return _value_text(values, limit)
    numbered = {}
    other = []
    for raw in dict.fromkeys(str(value) for value in values):
        match = re.fullmatch(r"(.*?)([0-9]+)", raw)
        if match and str(int(match[2])) == match[2]:
            numbered.setdefault(match[1], []).append(int(match[2]))
        else:
            other.append(raw)
    parts = []
    for prefix, numbers in numbered.items():
        runs = []
        for number in sorted(numbers):
            if runs and number == runs[-1][-1] + 1:
                runs[-1].append(number)
            else:
                runs.append([number])
        for run in runs:
            labels = ["%s%d" % (prefix, number) for number in run]
            compact = "%s–%s" % (labels[0], labels[-1])
            parts.extend([compact] if len(compact) < len(", ".join(labels)) else labels)
    return _value_text([*parts, *other], limit)


def _equivalent_metric(value, previous, settings):
    """Compare complete facts; never equate lists by a truncated screen prefix."""
    if isinstance(value, (list, tuple, dict)):
        try:
            return bool(value == previous)
        except (TypeError, ValueError):
            return False
    if isinstance(previous, (list, tuple)):
        return _value_text(value, settings.max_list_items) == _value_text(len(previous), settings.max_list_items)
    return _value_text(value, settings.max_list_items) == _value_text(previous, settings.max_list_items)


def _outcome_metric_key(label, records, settings):
    """Select only explicitly configured equivalents present for this file.

    In particular, a direct fallback record count must not be labelled as an
    index count when no index-count evidence exists.
    """
    choices = settings.outcome_metric_aliases.get(label, label)
    if isinstance(choices, str):
        return choices
    return next((key for key in choices if any(key in record.metrics for record in records)), choices[-1])


def _display_metric_key(key, value, records, settings):
    """Merge equivalent labels only when the corresponding evidence agrees."""
    target = _outcome_metric_key(key, records, settings)
    if target != key and any(target in record.metrics and _equivalent_metric(value, record.metrics[target], settings)
                             for record in records):
        return target
    return key


def _file_fields(group, settings, *, duplicate_name=False, file_number=None, update=False):
    records = group.records or group.all_records
    substantive = [record for record in records if set(record.checks) - _AVAILABILITY]
    named_records = [record for record in records if record.path is not None]
    preferred = max(enumerate(named_records or records), key=lambda item: (
        bool(set(item[1].checks) - _AVAILABILITY) or item[1].status != "passed",
        bool(item[1].metrics and set(item[1].metrics) != {"paths"}),
        bool(set(item[1].checks) - _AVAILABILITY), len(item[1].checks), item[0],
    ))[1]
    role = settings.role_labels.get(preferred.role, preferred.role)
    # Generic availability labels sometimes describe a failure condition.
    # They are audit labels, not safe headings for a successful file card.
    if not substantive and group.status == "passed":
        role = "File"
    state = group.status.upper()
    if group.status == "passed":
        state = "CHECKS PASSED" if substantive else (
            "AVAILABLE — availability only" if any(record.checks for record in records)
            else "OBSERVED — no check details"
        )
    title = "%s — %s" % (role, state)
    if update:
        title = "Additional findings for file %s — %s" % (file_number, state)
    elif file_number is not None:
        title = "%s · %s" % (file_number, title)
    fields = [screen_line(_STATUS_KIND[group.status], title, indent=6)]

    def add(kind, label, value, *, path_value=False):
        fields.append(screen_field(kind, label, value, indent=10,
                                   label_width=SUMMARY_CARD_LABEL_WIDTH, path_value=path_value))

    if not update:
        add("info", "File", group.path if duplicate_name else Path(group.path).name,
            path_value=duplicate_name)
    if group.companions:
        add("info", "Index", "; ".join(
            "%s: %s" % (Path(path).suffix.lstrip("."),
                         min((r.status for r in items), key=_STATUS_ORDER.__getitem__).upper())
            for path, items in group.companions.items()
        ))
    seen_metrics = set()
    for key, label in settings.metric_labels.items():
        values = [record.metrics[key] for record in records
                  if key in record.metrics and record.metrics[key] is not None]
        if not values:
            continue
        render = _contig_text if key in {"contigs", "chromosomes"} else _value_text
        for value in values:
            label = settings.metric_labels[_display_metric_key(key, value, records, settings)]
            text = render(value, settings.max_list_items)
            if (label, text) not in seen_metrics:
                add("count", label, text)
                seen_metrics.add((label, text))
    for record in records:
        for kind, label, value in record.metrics.get("reported_fields", ()):
            label = settings.metric_labels.get(_outcome_metric_key(label, records, settings), label)
            text = _value_text(value, settings.max_list_items)
            if (label, text) not in seen_metrics:
                add(kind, label, text)
                seen_metrics.add((label, text))
    checks = list(dict.fromkeys(check for record in reversed(substantive) for check in record.checks
                               if check not in _AVAILABILITY))
    if checks:
        text = "; ".join(checks[:settings.max_checks_per_file])
        if len(checks) > settings.max_checks_per_file:
            text += "; additional checks in audit"
        add("info", "Checks", text)
    messages = list(dict.fromkeys((record.status, record.message) for record in group.all_records
                                 if record.status != "passed" and record.message))
    for status, message in messages:
        if any(status == other_status and message != other and message in other
               for other_status, other in messages):
            continue
        add(_STATUS_KIND[status], status.capitalize(), message)
    return fields


def _validation_bundle_fields(bundle, settings, file_numbers):
    """Render only the caller-supplied summary of explicitly covered checks."""
    role = settings.role_labels[bundle.role]
    numbers = sorted(
        file_numbers[path] for path in bundle.paths if path in file_numbers
    )
    if numbers:
        label = "File" if len(numbers) == 1 else "Files"
        role = "%s %s · %s" % (
            label,
            _contig_text([str(number) for number in numbers], settings.max_list_items),
            role,
        )
    state = (
        "AVAILABLE — availability only"
        if bundle.availability_only
        else "CHECKS PASSED"
    )
    fields = [screen_line("success", "%s — %s" % (role, state), indent=6)]
    for kind, key, value, path_value in bundle.fields:
        render = _contig_text if key in {"contigs", "chromosomes"} else _value_text
        fields.append(screen_field(
            kind,
            settings.metric_labels[key],
            render(value, settings.max_list_items),
            indent=10,
            label_width=SUMMARY_CARD_LABEL_WIDTH,
            path_value=path_value,
        ))
    return fields


def render_file_validation(
    records,
    configuration,
    *,
    report_path=None,
    file_numbers=None,
    updated_paths=(),
    name_counts=None,
    omitted_count=0,
    validation_bundles=(),
    total_group_count=None,
    total_status_counts=None,
):
    """Return a shared, compact screen block; retain all details in the audit."""
    settings = configuration.logging.file_validation
    groups, requirements = group_file_validations(records)
    counts = (
        Counter(group.status for group in groups)
        if total_status_counts is None
        else Counter(total_status_counts)
    )
    group_count = len(groups) if total_group_count is None else total_group_count
    fields = ["", screen_line("analysis", "File validation", indent=2)]
    if group_count:
        fields.append(screen_field(
            "count", "Files (indexes grouped)",
            "%d — %s" % (group_count, "; ".join("%s: %d" % item for item in counts.items())),
            indent=6, label_width=SUMMARY_CARD_LABEL_WIDTH,
        ))
    names = name_counts if name_counts is not None else Counter(Path(group.path).name for group in groups)
    limit = settings.max_screen_files
    shown = [group for number, group in enumerate(groups)
             if limit is None or number < limit or group.status != "passed"]
    for group in shown:
        fields.extend(("", *_file_fields(
            group, settings, duplicate_name=names[Path(group.path).name] > 1,
            file_number=(file_numbers or {}).get(group.path),
            update=group.path in updated_paths,
        )))
    for bundle in validation_bundles:
        fields.extend(("", *_validation_bundle_fields(
            bundle, settings, file_numbers or {},
        )))
    shown_messages = {record.message for group in groups for record in group.all_records
                      if record.status != "passed" and record.message}
    for record in requirements:
        if record.status in {"failed", "blocked", "warning"} and record.message not in shown_messages:
            consumer = ", ".join(record.consumers)
            title = "%s%s" % (record.role, " (%s)" % consumer if consumer else "")
            fields.extend(("", screen_line(_STATUS_KIND[record.status], title, indent=6),
                           screen_field(_STATUS_KIND[record.status], record.status.capitalize(), record.message,
                                        indent=10, label_width=SUMMARY_CARD_LABEL_WIDTH)))
            shown_messages.add(record.message)
    deferred = sum(record.status == "deferred" for record in requirements)
    if deferred:
        fields.append(screen_field("decision", "Later checks", "%d deferred to consuming stages; see audit" % deferred,
                                   indent=6, label_width=SUMMARY_CARD_LABEL_WIDTH))
    omitted_count += len(groups) - len(shown)
    if omitted_count:
        fields.append(screen_field("info", "Additional files", "%d in the full audit" % omitted_count,
                                   indent=6, label_width=SUMMARY_CARD_LABEL_WIDTH))
    fields.extend(("", screen_line("info", "Only the listed checks are established; full details remain in the audit.", indent=6)))
    if report_path is not None:
        fields.append(screen_field("info", "Full validation audit", str(report_path), indent=6,
                                   label_width=SUMMARY_CARD_LABEL_WIDTH, path_value=True))
    return "\n".join((*fields, ""))


def _audit_destination(report_kind, records, path):
    """Check report ownership and input safety without creating a file."""
    destination = Path(path)
    if any(item.is_symlink() for item in (destination, *destination.parents)):
        raise ValueError("Validation report must not overwrite a symlink: %s" % path)
    destination = destination.resolve()
    for record in records:
        paths = (record.path,) if record.path is not None else record.metrics.get("paths", ())
        if any(Path(source).resolve() == destination for source in paths):
            raise ValueError("Validation report would overwrite a validated input: %s" % path)
    if destination.exists():
        expected = "report_kind: " + report_kind
        with destination.open("r", encoding="utf-8") as handle:
            marker = handle.readline(len(expected) + 2)
        if marker.rstrip("\r\n") != expected:
            raise ValueError("Refusing to overwrite an unrecognised validation audit file: %s" % path)
    return destination


def validation_audit_path(output_directory, pattern, **values):
    """Resolve the configured audit name without losing symlink evidence."""
    root = Path(output_directory).expanduser().resolve()
    destination = configured_output_path(root, str(pattern), **values)
    unresolved = root / str(pattern).format(**values)
    if any(path.is_symlink() for path in (unresolved, *unresolved.parents) if path != root):
        raise ValueError("Validation report must not use a symlink: %s" % unresolved)
    if destination == root:
        raise ValueError("Validation report must name a file, not the output root")
    return destination


def write_validation_audit(document, records, path):
    """Publish a small audit atomically without overwriting inputs or arbitrary files."""
    destination = _audit_destination(document["report_kind"], records, path)
    document = dict(document)
    document["report_version"] = 2
    display = getattr(current_validation_recorder(), "display", None)
    if display is not None:
        groups, _ = group_file_validations(records)
        for group in groups:
            display.file_numbers.setdefault(group.path, len(display.file_numbers) + 1)
    document["files"] = _combined_audit_records(records, getattr(display, "file_numbers", {}))
    return write_yaml_report(document, destination)


def _combined_audit_records(records, file_numbers=None):
    """Store each path once, with an indexed history of distinct check results.

    Metric values are interned per file. Check observations refer to their keys
    instead of repeating the path and facts. Conflicting values are retained,
    never silently overwritten or interpreted as equivalent checks.
    """
    by_path = {}
    requirements = []
    for record in records:
        if record.path is None:
            if record.to_dict() not in requirements:
                requirements.append(record.to_dict())
            continue
        by_path.setdefault(record.path, []).append(record)
    result = []
    for path, items in by_path.items():
        entry = dict(path=path, role=items[-1].role, roles=[], checks=[], status="passed",
                     metrics={}, message="", consumers=[], observations=[])
        if path in (file_numbers or {}):
            entry["file_number"] = file_numbers[path]
        messages = []
        for record in items:
            for key, values in (("roles", (record.role,)), ("checks", record.checks),
                                ("consumers", record.consumers)):
                entry[key].extend(value for value in values if value not in entry[key])
            if _STATUS_ORDER[record.status] < _STATUS_ORDER[entry["status"]]:
                entry["status"] = record.status
                entry["role"] = record.role
            if record.message and record.message not in messages:
                messages.append(record.message)
            metric_keys = []
            for key, value in record.metrics.items():
                identifier = key
                number = 1
                while identifier in entry["metrics"] and entry["metrics"][identifier] != value:
                    number += 1
                    identifier = "%s#%d" % (key, number)
                entry["metrics"][identifier] = value
                metric_keys.append(identifier)
            observation = dict(role=record.role, checks=list(record.checks), status=record.status,
                               metric_keys=metric_keys, consumers=list(record.consumers))
            if record.message:
                observation["message_index"] = messages.index(record.message)
            if observation not in entry["observations"]:
                entry["observations"].append(observation)
        entry["message"] = "\n".join(messages)
        entry["messages"] = messages
        result.append(entry)
    return [*result, *requirements]


def record_validation_fields(path, role, fields, *, checks=(), metrics=None):
    """Combine explicitly attributed module facts with the shared file record.

    Return the original fields when no reporter owns the invocation. Callers
    must supply file-specific facts, not analysis results or cross-file counts.
    No arbitrary terminal or native-tool text is filtered by this function.
    """
    if current_validation_recorder() is None:
        return fields
    values = dict(metrics or {})
    values["reported_fields"] = [list(item) for item in fields if len(item) == 3
                                and item[0] not in {"warning", "error", "loss"}]
    if values["reported_fields"] or metrics or checks:
        record_file_validation(path, role, checks=checks, metrics=values)
    for kind, label, value in (item for item in fields if len(item) == 3
                               and item[0] in {"warning", "error", "loss"}):
        record_file_validation(path, role, status="failed" if kind == "error" else "warning",
                               message="%s: %s" % (label, value))
    return []


def consolidate_validation_fields(fields, *, log_details=None):
    """Combine explicitly file-linked semantic fields, never arbitrary text.

    Exact paths (or unambiguous recorded basenames) establish attribution.
    Ambiguous and cross-file outcomes stay intact. Configured label equivalents
    remove only established facts; novel file facts join the same audit record.
    This function does not read files, infer check success, or alter validation.
    """
    session = current_validation_recorder()
    display = getattr(session, "display", None)
    if not fields or display is None:
        return fields
    records = session.records
    paths = {record.path for record in records if record.path is not None}
    groups = [[]]
    for item in fields:
        if len(item) == 2 and groups[-1]:
            groups.append([])
        groups[-1].append(item)
    retained = []
    attributed_paths = set()
    settings = display.configuration.logging.file_validation
    for group in groups:
        candidates = set()
        path_fields = []
        references = {}
        for item in group:
            if len(item) != 3 or not isinstance(item[2], (str, Path)):
                continue
            value = str(item[2])
            matches = {path for path in paths if value == path}
            if not matches:
                matches = {path for path in paths if value == Path(path).name}
            if len(matches) == 1:
                candidates.update(matches)
                path_fields.append(item)
                references[value] = next(iter(matches))
        if len(candidates) > 1:
            # Keep cross-file findings at their consuming stage. Refer to the
            # existing file sections instead of repeating their inventories.
            retained.extend(
                (item[0], item[1], display.file_reference(references[str(item[2])]))
                if item in path_fields else item for item in group
            )
            continue
        if len(candidates) != 1:
            retained.extend(group)
            continue
        path = candidates.pop()
        attributed_paths.add(path)
        facts = [record for record in records if record.path == path]
        novel = []
        for item in group:
            if len(item) == 2 or item in path_fields:
                continue
            kind, label, value = item
            key = _outcome_metric_key(label, facts, settings)
            established = [record.metrics[key] for record in facts if key in record.metrics] if key else []
            equivalent = any(_equivalent_metric(value, previous, settings) for previous in established)
            if kind in {"warning", "error", "loss"} or not equivalent:
                novel.append(item)
        record_validation_fields(path, facts[-1].role, novel)
    if log_details is not None and not retained and len(attributed_paths) == 1:
        path = next(iter(attributed_paths))
        # Only the copy intended for terminal/file-log presentation is moved.
        # StepContext.extra and its downstream scientific consumers retain the
        # original structured values unchanged.
        facts = [record for record in session.records if record.path == path]
        metrics = {_outcome_metric_key(key, facts, settings): value
                   for key, value in log_details.items()
                   if not isinstance(value, (str, Path)) or str(value) not in {path, Path(path).name}}
        metrics = {key: value for key, value in metrics.items() if not any(
            key in record.metrics and _equivalent_metric(value, record.metrics[key], settings)
            for record in session.records if record.path == path
        )}
        try:
            if metrics:
                record_file_validation(path, "Stage validation details", metrics=metrics)
        except (TypeError, ValueError):
            # Rich analysis objects and non-finite scientific diagnostics are
            # not compact validation evidence. Keep their original log output;
            # presentation must never make a valid analysis fail.
            pass
        else:
            log_details.clear()
    return retained


def consolidate_validation_log(marker, subject, values):
    """Combine explicitly classified validation events, not free-form tool logs."""
    session = current_validation_recorder()
    display = getattr(session, "display", None)
    if marker != "VALIDATE" or display is None:
        return values
    paths = {record.path for record in session.records if record.path is not None}
    matched = {str(value) for value in values.values()
               if isinstance(value, (str, Path)) and str(value) in paths}
    if len(matched) != 1:
        return values
    path = matched.pop()
    status = str(values.get("status", "")).lower()
    if status not in _STATUS_ORDER:
        return values
    metrics = {key: value for key, value in values.items()
               if key not in {"status", "message"}
               and not (isinstance(value, (str, Path)) and str(value) == path)}
    metrics = {key: value for key, value in metrics.items() if not any(
        key in record.metrics and _equivalent_metric(
            value, record.metrics[key], display.configuration.logging.file_validation,
        )
        for record in session.records if record.path == path
    )}
    try:
        record_file_validation(path, str(subject), status=status, metrics={"logged_validation": metrics},
                               message=str(values.get("message") or ""))
    except (TypeError, ValueError):
        return values
    if not display._new((path, "log", subject), values):
        return None
    return {"validation_file": display.file_reference(path), "status": values["status"]}


class FileValidationDisplay:
    """Invocation-local, incremental presentation of established file facts."""

    def __init__(self, session, configuration):
        self.session = session
        self.configuration = configuration
        self.owner_pid = os.getpid()
        self.file_numbers = {}
        self.displayed_paths = set()
        self.seen = set()
        self.report_path_displayed = False
        self.validation_bundles = {}
        session.display = self

    def _new(self, key, value):
        fact = (key, json.dumps(value, sort_keys=True, default=str))
        if fact in self.seen:
            return False
        self.seen.add(fact)
        return True

    def file_reference(self, path):
        """Number an observed file without rescanning it or guessing its role."""
        groups, _ = group_file_validations(self.session.records)
        for group in groups:
            number = self.file_numbers.setdefault(group.path, len(self.file_numbers) + 1)
            if path == group.path or path in group.companions:
                return "File %s" % number
        raise ValueError("Cannot reference an unrecorded validation file")

    def register_validation_bundle(
        self,
        paths,
        role,
        fields,
        *,
        covered_checks=(),
        covered_metric_keys=(),
        availability_only,
    ):
        """Register compact screen metadata without changing validation evidence."""
        resolved_paths = tuple(dict.fromkeys(
            str(Path(path).expanduser().resolve()) for path in paths
        ))
        if not resolved_paths:
            raise ValueError("A file-validation bundle requires at least one path")
        settings = self.configuration.logging.file_validation
        if role not in settings.role_labels:
            raise ValueError("Unknown configured file-validation role: %s" % role)
        normalized_fields = []
        for item in fields:
            if len(item) not in {3, 4}:
                raise ValueError(
                    "File-validation bundle fields require three or four values"
                )
            kind, key, value = item[:3]
            if kind not in self.configuration.logging.terminal_style.styles:
                raise ValueError("Unknown configured terminal role: %s" % kind)
            if key not in settings.metric_labels:
                raise ValueError("Unknown configured file-validation metric: %s" % key)
            normalized_fields.append((
                str(kind), str(key), value,
                bool(item[3]) if len(item) == 4 else False,
            ))
        bundle = FileValidationBundle(
            str(role),
            resolved_paths,
            tuple(normalized_fields),
            frozenset(str(value) for value in covered_checks),
            frozenset(str(value) for value in covered_metric_keys),
            bool(availability_only),
        )
        key = (bundle.role, bundle.paths)
        previous = self.validation_bundles.get(key)
        if previous is not None and previous != bundle:
            raise ValueError(
                "Conflicting summaries registered for file-validation bundle %s"
                % role
            )
        self.validation_bundles[key] = bundle

    def register_availability_bundle(
        self, paths, role, fields, *, covered_metric_keys=(),
    ):
        """Register a bundle whose successful claim is availability only."""
        self.register_validation_bundle(
            paths,
            role,
            fields,
            covered_metric_keys=covered_metric_keys,
            availability_only=True,
        )

    @staticmethod
    def _bundle_member_is_covered(group, bundle):
        """Require every hidden fact to be represented by the declared bundle."""
        if group.status != "passed" or not group.all_records:
            return False
        has_check = False
        has_substantive_check = False
        for record in group.all_records:
            checks = set(record.checks)
            has_check = has_check or bool(checks)
            has_substantive_check = has_substantive_check or bool(
                checks - _AVAILABILITY
            )
            if checks - _AVAILABILITY - bundle.covered_checks:
                return False
            if set(record.metrics) - bundle.covered_metric_keys:
                return False
        if bundle.availability_only:
            return has_check and not has_substantive_check
        return has_substantive_check

    def _mark_bundle_group_displayed(self, group):
        """Prevent a summarized success from being repeated on a later flush."""
        settings = self.configuration.logging.file_validation
        for record in group.all_records:
            identity = record.path or group.path
            for check in record.checks:
                self._new((identity, "check", record.status), check)
            for key, value in record.metrics.items():
                if key == "reported_fields":
                    for kind, label, fact in value:
                        label = settings.metric_labels.get(
                            _outcome_metric_key(label, group.records, settings),
                            label,
                        )
                        identity_value = (
                            fact
                            if isinstance(fact, (list, tuple, dict))
                            else _value_text(fact, settings.max_list_items)
                        )
                        self._new((identity, "display", label), identity_value)
                elif key in settings.metric_labels:
                    label = settings.metric_labels[
                        _display_metric_key(
                            key, value, group.records, settings,
                        )
                    ]
                    identity_value = (
                        value
                        if isinstance(value, (list, tuple, dict))
                        else _value_text(value, settings.max_list_items)
                    )
                    self._new((identity, "display", label), identity_value)
            if record.message:
                self._new((identity, "message"), record.message)
            self._new((identity, "status"), record.status)
        self.displayed_paths.add(group.path)

    def flush(self, *, report_path=None):
        if os.getpid() != self.owner_pid:
            return
        groups, requirements = group_file_validations(self.session.records)
        for group in groups:
            self.file_numbers.setdefault(group.path, len(self.file_numbers) + 1)
        pending = []
        updated = set(self.displayed_paths)
        settings = self.configuration.logging.file_validation
        omitted_count = 0
        groups_by_path = {group.path: group for group in groups}
        pending_bundles = []
        bundled_paths = set()
        for key, bundle in self.validation_bundles.items():
            members = [groups_by_path.get(path) for path in bundle.paths]
            if (
                not (set(bundle.paths) & bundled_paths)
                and all(
                    member is not None
                    and self._bundle_member_is_covered(member, bundle)
                    for member in members
                )
                and self._new(("validation_bundle", key), bundle.fields)
            ):
                pending_bundles.append(bundle)
                bundled_paths.update(bundle.paths)
                for member in members:
                    self._mark_bundle_group_displayed(member)

        presentation_number = 0
        for group in groups:
            if group.path in bundled_paths:
                continue
            if (settings.max_screen_files is not None and presentation_number >= settings.max_screen_files
                    and group.status == "passed"):
                omitted_count += 1
                presentation_number += 1
                continue
            presentation_number += 1
            group_start = len(pending)
            for record in group.all_records:
                identity = record.path or group.path
                checks = tuple(check for check in record.checks
                               if self._new((identity, "check", record.status), check))
                metrics = {}
                for key, value in record.metrics.items():
                    if key == "reported_fields":
                        fields = []
                        for item in value:
                            kind, label, fact = item
                            label = settings.metric_labels.get(_outcome_metric_key(label, group.records, settings), label)
                            identity_value = fact if isinstance(fact, (list, tuple, dict)) else _value_text(fact, settings.max_list_items)
                            if self._new((identity, "display", label), identity_value):
                                fields.append(item)
                        if fields:
                            metrics[key] = fields
                    elif key in settings.metric_labels and self._new(
                        (identity, "display", settings.metric_labels[
                            _display_metric_key(key, value, group.records, settings)
                        ]),
                        value if isinstance(value, (list, tuple, dict)) else _value_text(value, settings.max_list_items),
                    ):
                        metrics[key] = value
                message = record.message if record.message and self._new((identity, "message"), record.message) else ""
                first_status = self._new((identity, "status"), record.status)
                if checks or metrics or message or first_status:
                    pending.append(replace(record, checks=checks, metrics=metrics, message=message))
            if len(pending) > group_start:
                if group.records and not any(record.path == group.path for record in pending[group_start:]):
                    pending.append(replace(group.records[0], checks=(), metrics={}, message=""))
                self.displayed_paths.add(group.path)
        for record in requirements:
            value = record.to_dict()
            value.pop("consumers", None)
            if self._new("requirement", value):
                pending.append(record)
        if pending or pending_bundles:
            print_screen_block(render_file_validation(
                pending, self.configuration,
                report_path=report_path if not self.report_path_displayed else None,
                file_numbers=self.file_numbers, updated_paths=updated,
                name_counts=Counter(Path(group.path).name for group in groups),
                omitted_count=omitted_count,
                validation_bundles=pending_bundles,
                total_group_count=len(groups) if pending_bundles else None,
                total_status_counts=(
                    Counter(group.status for group in groups)
                    if pending_bundles else None
                ),
            ))
            if report_path is not None:
                self.report_path_displayed = True


@contextmanager
def file_validation_display(session, configuration, *, display=None):
    """Use the same stage-boundary display in direct and pipeline execution."""
    display = display or FileValidationDisplay(session, configuration)
    token = _DIRECT_DISPLAY.set(display.flush)
    try:
        yield display
    finally:
        _DIRECT_DISPLAY.reset(token)


def flush_file_validation_display():
    """Flush invocation evidence at progress boundaries, without new checks."""
    callback = _DIRECT_DISPLAY.get()
    if callback is not None:
        callback()


def register_file_availability_bundle(
    paths, role, fields, *, covered_metric_keys=(),
):
    """Request compact display of an explicit bundle in the active invocation.

    This is presentation metadata only. Every exact file record remains in the
    canonical audit, and the display uses the bundle only when all members have
    passed availability checks with no stronger claims. Any per-file metrics
    must be explicitly named and represented by the supplied summary fields.
    """
    session = current_validation_recorder()
    display = getattr(session, "display", None)
    if display is not None:
        display.register_availability_bundle(
            paths,
            role,
            fields,
            covered_metric_keys=covered_metric_keys,
        )


def register_file_validation_bundle(
    paths,
    role,
    fields,
    *,
    covered_checks,
    covered_metric_keys=(),
):
    """Request compact display of an explicitly summarized validated bundle.

    Every member and check remains in the audit. The bundle is rendered only
    when every member passed and every recorded substantive check and metric is
    declared here; otherwise ordinary per-file cards remain visible.
    """
    session = current_validation_recorder()
    display = getattr(session, "display", None)
    if display is not None:
        display.register_validation_bundle(
            paths,
            role,
            fields,
            covered_checks=covered_checks,
            covered_metric_keys=covered_metric_keys,
            availability_only=False,
        )


def run_with_file_validation(operation, *, command, configuration, output_directory):
    """Reuse exact direct checks without enabling pipeline-only contracts.

    The report is written after the operation. An orchestrated checkpoint wraps
    this function; native module checkpoints never own an in-progress audit.
    Commands that perform no recorded checks (including a reused checkpoint)
    do not invent a fresh validation pass or overwrite a previous audit.
    """
    if current_validation_session() is not None:
        return operation()
    validate_filename_component(command, "direct command")
    destination = validation_audit_path(output_directory, configuration.logging.file_validation.direct_report_file,
                                        command=command)
    _audit_destination(_DIRECT_REPORT_KIND, (), destination)
    with (InputValidationSession(pipeline=False) as session, session.scope(command),
          file_validation_display(session, configuration) as display):
        failure = None
        try:
            result = operation()
            if isinstance(result, int) and not isinstance(result, bool) and result != 0:
                failure = RuntimeError("Command returned exit status %d" % result)
            return result
        except BaseException as exc:
            if not isinstance(exc, SystemExit) or exc.code not in (None, 0):
                failure = exc
            raise
        finally:
            try:
                if session.records:
                    if failure is not None:
                        record_file_validation(None, "Command execution", status="failed", message=str(failure))
                    display.flush()
                    write_validation_audit({
                        "report_kind": _DIRECT_REPORT_KIND,
                        "phase": "direct_execution", "command": command,
                        "status": "failed" if failure is not None else "passed",
                        "error": None if failure is None else "%s: %s" % (type(failure).__name__, failure),
                        "interpretation": "Only observed checks are reported; no pipeline-only requirements or additional scans were imposed.",
                    }, session.records, destination)
                    print_screen_block(screen_field("info", "Full validation audit", str(destination), indent=6,
                                                    label_width=SUMMARY_CARD_LABEL_WIDTH, path_value=True))
            except Exception as report_error:
                if failure is None:
                    raise
                print_screen_block(screen_line("error", "Validation audit could not be saved: %s" % report_error))


__all__ = ["group_file_validations", "render_file_validation", "write_validation_audit",
           "validation_audit_path", "flush_file_validation_display", "run_with_file_validation",
           "file_validation_display", "FileValidationDisplay", "record_validation_fields",
           "consolidate_validation_fields", "consolidate_validation_log",
           "register_file_availability_bundle", "register_file_validation_bundle"]
