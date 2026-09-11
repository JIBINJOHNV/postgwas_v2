"""Canonical structured logging for every PostGWAS module and pipeline scope.

Every PostGWAS component can use this module so the audit trail has one
stable, searchable vocabulary:

  1. which chromosome, on every line          -> the scope column
  2. which step and which function            -> step()
  3. variant count on entry                   -> step(rows_in=...)
  4. every policy value actually consulted    -> PARAM records
  5. every decision and effective QC action    -> OBSERVED / ACTION / RESULT
     with before and after counts                 records
  6. outcome, STATUS or FAILED with the        -> step() logs both, and re-raises
     plain-English cause then the traceback

Line format (section 3.2):

    [YYYY-MM-DD HH:MM:SS] chrN   MARKER  message

The file log is compact and machine-searchable.  Long policy explanations live
in the canonical YAML and ``postgwas config export --style full``; repeating
them for every step and chromosome would obscure what the analysis did.
Zero-effect QC checks are retained as one-line PASS records.  Explanations and
actions are emitted only when a check affected data or raised a warning.

Workers never write to stdout.  Each buffers its own screen text; the parent
prints one complete block per chromosome with print_chromosome_summary(), from its
single-threaded loop, so interleaving is structurally impossible.

The log file is opened in APPEND mode and flushed on every record, so a SIGKILL
still leaves a complete log and every retry of a chromosome accumulates in one
file.

The same formatter is used for run, dataset, chromosome, and concordance logs.
"""

import datetime
import os
import re
import sys
import textwrap
import threading
import time
import traceback
from contextlib import contextmanager
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from postgwas.core.ui.screen import (
    normalize_stage_outcome_fields,
    screen_field,
    screen_line,
)
from postgwas.core.values import format_count as _count

__all__ = [
    "PipelineLogger",
    "StepContext",
    "write_log_record",
    "print_chromosome_summary",
    "MARKERS",
    "LEVELS",
]


def write_log_record(
    path,
    level,
    message,
    *,
    sample_id="postgwas",
    scope="run",
    file_level="DEBUG",
    screen_level="ERROR",
):
    """Write and flush one canonical record, including during preflight failure."""
    logger = PipelineLogger(
        sample_id=sample_id,
        scope=scope,
        log_dir=os.path.dirname(os.path.abspath(os.fspath(path))),
        level=file_level,
        screen_level=screen_level,
        log_path=os.fspath(path),
    )
    try:
        normalized = str(level).upper()
        if normalized == "DEBUG":
            logger.debug(message)
        elif normalized == "WARNING":
            logger.warning(message)
        elif normalized in ("ERROR", "CRITICAL"):
            logger.error(message)
        else:
            logger.info(message)
    finally:
        logger.close()


# =============================================================================
# Constants
# =============================================================================

MARKERS = (
    "STEP", "INPUT", "PARAM", "OBSERVED", "DECIDE", "ACTION", "RESULT",
    "OUTPUT", "STATUS", "PASS", "DONE", "FAILED", "WARNING", "SKIP", "",
)

LEVELS = {"DEBUG": 10, "INFO": 20, "WARNING": 30, "ERROR": 40}

_MARKER_LEVEL = {
    "STEP": "INFO",
    "INPUT": "INFO",
    "PARAM": "INFO",
    "OBSERVED": "INFO",
    "ACTION": "INFO",
    "RESULT": "INFO",
    "OUTPUT": "INFO",
    "STATUS": "INFO",
    "PASS": "INFO",
    "DONE": "INFO",
    "SKIP": "INFO",
    "DECIDE": "INFO",
    "WARNING": "WARNING",
    "FAILED": "ERROR",
    "": "INFO",
}

_TIMESTAMP_WIDTH = 22   # "[2026-08-03 17:59:01] "
_SCOPE_WIDTH = 7        # "chr7   "
_MARKER_WIDTH = 9       # "OBSERVED " (one separator after the longest marker)
_LINE_WIDTH = 128       # where messages and help text are wrapped

_PREFIX_WIDTH = _TIMESTAMP_WIDTH + _SCOPE_WIDTH + _MARKER_WIDTH


# =============================================================================
# Small formatting helpers
# =============================================================================


def _elapsed(seconds):
    # type: (Optional[float]) -> str
    if seconds is None:
        return "unknown"
    seconds = float(seconds)
    if seconds < 60:
        return "%.1f s" % seconds
    if seconds < 3600:
        return "%d m %02d s" % (int(seconds // 60), int(seconds % 60))
    return "%d h %02d m" % (int(seconds // 3600), int((seconds % 3600) // 60))


def _format_value(value):
    # type: (Any) -> str
    """Render a policy value.  Matches policies.format_value; kept local so that
    this module has no imports of its own beyond the standard library."""
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if value != value:
            return "nan"
        if value == 0.0:
            return "0.0"
        if abs(value) < 1e-4 or abs(value) >= 1e6:
            return repr(value)
        text = ("%.10f" % value).rstrip("0")
        if text.endswith("."):
            text += "0"
        return text
    if isinstance(value, (list, tuple)):
        return "[" + ", ".join(_format_value(item) for item in value) + "]"
    if isinstance(value, dict):
        return "{" + ", ".join(
            "%s -> %s" % (key, _format_value(item)) for key, item in value.items()
        ) + "}"
    return str(value)


def _field_value(value):
    # type: (Any) -> str
    """Render one structured field without making paths hard to copy."""
    text = _format_value(value)
    if text == "" or any(character.isspace() for character in text):
        return repr(text)
    return text


def _fields(**values):
    return " ".join(
        "%s=%s" % (name, _field_value(value))
        for name, value in values.items() if value is not None
    )


def _wrap_lines(lines, width):
    # type: (Sequence[str], int) -> List[str]
    """Wrap at word boundaries, keeping each line's own leading spaces so that
    aligned blocks (the settings table) stay aligned.  Long single words such as
    file paths are never broken."""
    out = []  # type: List[str]
    for line in lines:
        if len(line) <= width:
            out.append(line)
            continue
        stripped = line.lstrip(" ")
        lead = " " * (len(line) - len(stripped))
        pieces = textwrap.wrap(
            stripped,
            max(20, width - len(lead)),
            break_long_words=False,
            break_on_hyphens=False,
        )
        if not pieces:
            out.append(line)
            continue
        out.extend(lead + piece for piece in pieces)
    return out


def _step_number(number, total):
    # type: (Any, Any) -> str
    """'4', 16 -> '04/16', matching the plan's step labels."""
    try:
        width = len(str(int(total)))
        return "%0*d/%d" % (width, int(number), int(total))
    except (TypeError, ValueError):
        return "%s/%s" % (number, total)


def _scope_label(scope):
    # type: (str) -> str
    """'7' -> 'chr7', 'MT' -> 'chrMT', 'dataset' -> 'DS'."""
    text = str(scope).strip()
    if text.lower() in ("dataset", "ds", ""):
        return "DS"
    if text.lower() == "run":
        return "RUN"
    if text.lower() in ("validation", "concordance"):
        return "VALID"
    if text.lower().startswith("chr"):
        return "chr" + text[3:]
    return "chr" + text


def _is_dataset(scope):
    # type: (str) -> bool
    return str(scope).strip().lower() in ("dataset", "ds", "")


# =============================================================================
# Plain-English explanations for failures
# =============================================================================

_ERROR_EXPLANATIONS = [
    (
        "FileNotFoundError",
        "A file this step needs is not there.",
        "Check the paths in the config, and check that an earlier step actually "
        "produced the file this one is reading.",
    ),
    (
        "PermissionError",
        "The pipeline is not allowed to read or write one of the files this step uses.",
        "Check the permissions on the output folder and on the input files.",
    ),
    (
        "IsADirectoryError",
        "A path that should name a file names a folder instead.",
        "Check the file paths in the config.",
    ),
    (
        "MemoryError",
        "The machine ran out of memory while this step was working.",
        "Give the job more memory, or run fewer chromosomes at once.",
    ),
    (
        "ColumnNotFoundError",
        "This step asked for a column that is not in the data.",
        "Either the config names a column the file does not have, or an earlier "
        "step did not create the column this one expects.",
    ),
    (
        "DuplicateError",
        "Two columns ended up with the same name, so the data could not be read.",
        "Two config keys probably point at the same column in the input file, or "
        "the file itself has the same column name twice.",
    ),
    (
        "SchemaError",
        "A column has a different type from the one this step expected.",
        "A numeric column probably arrived as text because at least one value in "
        "it could not be read as a number.",
    ),
    (
        "SchemaFieldNotFoundError",
        "This step asked for a field that is not in the data.",
        "An earlier step did not create the column this one expects.",
    ),
    (
        "InvalidOperationError",
        "A calculation could not be carried out on the data as it stands.",
        "This is usually a column that is text where a number was expected.",
    ),
    (
        "ComputeError",
        "A calculation failed part-way through the data.",
        "This is usually one unreadable value in an otherwise numeric column.",
    ),
    (
        "ZeroDivisionError",
        "A calculation divided by zero.",
        "A frequency of exactly 0 or 1, or a sample size of 0, is the usual cause.",
    ),
    (
        "KeyError",
        "Something this step looked up was not there.",
        "Usually a missing config key or a missing column.",
    ),
    (
        "TimeoutError",
        "An external tool took longer than it was allowed and was stopped.",
        "Re-run the chromosome; if it happens again the input for this step is "
        "probably much larger than expected.",
    ),
    (
        "CalledProcessError",
        "An external tool (bcftools, tabix or similar) reported an error.",
        "The tool's own message is in the developer section below.",
    ),
]


def _explain_exception(exc):
    # type: (BaseException) -> Tuple[str, str]
    """(what went wrong, what to do) in plain English."""
    hint = getattr(exc, "plain_english", None)
    todo = getattr(exc, "what_to_do", None)
    if hint:
        return str(hint), str(todo or "See the developer section below.")

    names = [type(exc).__name__] + [base.__name__ for base in type(exc).__mro__]
    for name, what, action in _ERROR_EXPLANATIONS:
        if name in names:
            return what, action
    return (
        "This step stopped with an error it did not expect.",
        "The Python error and the traceback below say where; they are for developers.",
    )


# =============================================================================
# Step context
# =============================================================================


class StepContext(object):
    """The small mutable object step() yields.

    Set `rows_out` (and optionally `removed`) so the DONE line can report what
    the step did. Set one concise ``outcome`` for the canonical RESULT record;
    optional semantic fields render as an aligned terminal block. When
    `rows_out` is left alone it is taken to be unchanged. `extra` is carried
    into summary() for the manifest.
    """

    __slots__ = (
        "number", "total", "title", "func_name",
        "rows_in", "rows_out", "removed", "extra",
        "outcome_fields", "failure_hint", "_logger", "started", "logged_outcome",
    )

    def __init__(self, logger, number, total, title, func_name, rows_in):
        self.number = number
        self.total = total
        self.title = title
        self.func_name = func_name
        self.rows_in = rows_in
        self.rows_out = None    # type: Optional[int]
        self.removed = None     # type: Optional[int]
        self.extra = {}         # type: Dict[str, Any]
        self.outcome_fields = None
        self.logged_outcome = None
        self.failure_hint = None  # type: Optional[str]
        self.started = time.time()
        self._logger = logger

    # convenience passthroughs, so a step never has to reach past its context
    def info(self, message, indent=0):
        self._logger.info(message, indent=indent)

    def warn(self, message, indent=0):
        self._logger.warn(message, indent=indent)

    def error(self, message, indent=0):
        self._logger.error(message, indent=indent)

    def skip(self, message, indent=0):
        self._logger.skip(message, indent=indent)

    def decide(self, what, evidence, decision):
        self._logger.decide(what, evidence, decision)

    def qc(self, check_name, plain_english, before, after, **kwargs):
        self._logger.qc(check_name, plain_english, before, after, **kwargs)

    def input(self, subject, **values):
        self._logger.record("INPUT", subject, **values)

    def observed(self, subject, **values):
        self._logger.record("OBSERVED", subject, **values)

    def action(self, subject, **values):
        self._logger.record("ACTION", subject, **values)

    def output(self, subject, **values):
        self._logger.record("OUTPUT", subject, **values)

    def set_rows(self, rows_out, removed=None):
        self.rows_out = rows_out
        if removed is not None:
            self.removed = removed

    def outcome(self, message, fields=None, **values):
        from postgwas.core.validation_reporting import consolidate_validation_fields

        text = str(message).strip()
        if not text:
            raise ValueError("A completed-stage outcome must not be empty")
        reserved = {"message", "step"}.intersection(values)
        if reserved:
            raise ValueError(
                "Stage outcome details use reserved names: %s"
                % ", ".join(sorted(reserved))
            )
        self.extra["outcome"] = {"message": text, **values}
        log_values = dict(values)
        self.outcome_fields = normalize_stage_outcome_fields(consolidate_validation_fields(
            normalize_stage_outcome_fields(fields), log_details=log_values,
        ))
        self.logged_outcome = {"message": text, **log_values}

    @property
    def label(self):
        return "%s  %s" % (_step_number(self.number, self.total), self.title)

    def __repr__(self):
        return "<StepContext %s>" % self.label


# =============================================================================
# The logger
# =============================================================================


class PipelineLogger(object):
    """Structured logger shared by module, dataset, and chromosome scopes.

    log file:  {log_dir}/{sample_id}_chr{scope}.log   for a chromosome
               {log_dir}/{sample_id}_dataset.log      for scope == 'dataset'

    ``log_path`` may select an exact run-level or validation log. Otherwise the
    standard dataset/chromosome filename is derived from ``sample_id`` and
    ``scope``. Files open in append mode and flush after every record.
    """

    def __init__(
        self,
        sample_id,
        scope,
        log_dir,
        policies=None,
        level="INFO",
        screen_level="WARNING",
        log_path=None,
        stage_progress=None,
    ):
        # type: (str, str, str, Any, Optional[str], Optional[str], Optional[str], Any) -> None
        self.sample_id = str(sample_id)
        self.scope = str(scope)
        self.policies = policies
        self.label = _scope_label(scope)
        self.is_dataset = _is_dataset(scope)

        self.level = self._resolve_level(level, "logging.level", "INFO")
        self.screen_level = self._resolve_level(screen_level, "logging.screen_level", "WARNING")

        if log_path is not None:
            self.log_path = os.path.abspath(os.fspath(log_path))
            self.log_dir = os.path.dirname(self.log_path) or os.curdir
        else:
            self.log_dir = str(log_dir)
            if self.is_dataset:
                filename = "%s_dataset.log" % self.sample_id
            else:
                filename = "%s_chr%s.log" % (self.sample_id, self.scope)
            self.log_path = os.path.join(self.log_dir, filename)
        os.makedirs(self.log_dir, exist_ok=True)

        self._lock = threading.Lock()
        self._fh = None
        self._closed = False
        self._open()

        self._screen = []          # type: List[str]
        self._steps = []           # type: List[Dict[str, Any]]
        self._qc_actions = []      # type: List[Dict[str, Any]]
        self._decisions = []       # type: List[Dict[str, Any]]
        self._current_step = None  # type: Optional[StepContext]
        self.warning_count = 0
        self.error_count = 0
        self.started = time.time()
        self._stage_progress = stage_progress

    # -- plumbing ---------------------------------------------------------
    def _resolve_level(self, value, policy_key, fallback):
        if value is None and self.policies is not None:
            try:
                value = self.policies.get(policy_key)
            except Exception:
                value = None
        if value is None:
            value = fallback
        text = str(value).strip().upper()
        if text not in LEVELS:
            text = fallback
        return text

    def _open(self):
        if self._fh is None or self._fh.closed:
            # append mode: retries accumulate in one file
            self._fh = open(self.log_path, "a", encoding="utf-8")
            self._closed = False

    def __getstate__(self):
        state = dict(self.__dict__)
        state["_fh"] = None
        state["_lock"] = None
        # Live terminal state belongs only to the parent process.
        state["_stage_progress"] = None
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self._lock = threading.Lock()
        self._fh = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()
        return False

    def _progress_event(self, method, *values):
        observer = self._stage_progress
        if observer is None:
            return
        try:
            getattr(observer, method)(*values)
        except Exception as exc:
            self._stage_progress = None
            try:
                observer.close()
            except Exception:
                pass
            self.record(
                "WARNING",
                "terminal_progress_disabled",
                error="%s: %s" % (type(exc).__name__, exc),
            )

    def set_stage_progress(self, stage_progress):
        """Replace the live stage observer at a validated phase boundary.

        A logger may contain several independently numbered stage groups. The
        previous observer is closed before the next one is installed so a
        completed live region cannot remain active or receive events from a
        later group whose numbering starts again at one.
        """
        previous = self._stage_progress
        if previous is stage_progress:
            return
        self._stage_progress = None
        if previous is not None:
            try:
                previous.close()
            except Exception as exc:
                self.record(
                    "WARNING",
                    "terminal_progress_close_failed",
                    error="%s: %s" % (type(exc).__name__, exc),
                )
        self._stage_progress = stage_progress

    # -- the one and only line formatter ----------------------------------
    def log(self, marker, message, indent=0, level=None, wrap=True, screen=True):
        # type: (str, Any, int, Optional[str], bool, bool) -> None
        """Format and record one record.  Every other method calls this one.

        A message containing newlines becomes several lines: the marker appears
        on the first, the rest are indented detail lines underneath it.  Long
        lines are wrapped at word boundaries and continue as detail lines, with
        any leading spaces preserved so aligned blocks stay aligned.  Pass
        wrap=False for text that must not be reflowed, such as a traceback.
        """
        marker = (marker or "").upper()
        if marker not in _MARKER_LEVEL:
            marker = ""
        if level is None:
            level = _MARKER_LEVEL[marker]
        level = level.upper()

        if marker == "WARNING":
            self.warning_count += 1
        elif marker == "FAILED":
            self.error_count += 1

        stamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        pad = "  " * max(0, int(indent))
        lines = str(message).split("\n")
        if wrap:
            lines = _wrap_lines(lines, max(40, _LINE_WIDTH - _PREFIX_WIDTH - len(pad)))

        rendered = []
        for position, text in enumerate(lines):
            shown_marker = marker if position == 0 else ""
            rendered.append(
                "[%s] %-*s%-*s%s%s"
                % (
                    stamp,
                    _SCOPE_WIDTH, self.label,
                    _MARKER_WIDTH, shown_marker,
                    pad,
                    text,
                )
            )
        block = "\n".join(rendered)

        if LEVELS[level] >= LEVELS[self.level]:
            self._write(block)
        if screen and LEVELS[level] >= LEVELS[self.screen_level]:
            self._screen.extend(rendered)

    def _write(self, block):
        lock = self._lock or threading.Lock()
        with lock:
            self._open()
            self._fh.write(block + "\n")
            self._fh.flush()   # a SIGKILL still leaves a complete log

    # -- everyday records --------------------------------------------------
    def info(self, message, indent=0):
        self.log("", message, indent=indent, level="INFO")

    def debug(self, message, indent=0):
        self.log("", message, indent=indent, level="DEBUG")

    def warn(self, message, indent=0):
        self.log("WARNING", message, indent=indent)

    warning = warn

    def error(self, message, indent=0):
        self.log("FAILED", message, indent=indent)

    def skip(self, message, indent=0):
        self.log("SKIP", message, indent=indent)

    def blank(self):
        self.log("", "")

    def record(self, marker, subject, **values):
        """Write one compact structured audit record."""
        from postgwas.core.validation_reporting import consolidate_validation_log

        values = consolidate_validation_log(marker, subject, values)
        if values is None:
            return
        detail = _fields(**values)
        message = str(subject) if not detail else "%s %s" % (subject, detail)
        self.log(marker, message)

    # -- settings block ----------------------------------------------------
    def settings(self, policy_keys, indent=0):
        # type: (Iterable[str], int) -> None
        """Record resolved values only; canonical YAML owns the documentation."""
        keys = list(policy_keys or [])
        if not keys:
            return
        if self.policies is None:
            self.warn(
                "PARAM records unavailable: the logger received no policy registry. "
                "This is an internal wiring error, not a data problem."
            )
            return

        for key in keys:
            value = self.policies.get(key)
            self.log("PARAM", "%s=%s" % (key, _field_value(value)), indent=indent)

    # -- decisions ---------------------------------------------------------
    def decide(self, what, evidence, decision):
        # type: (str, Any, str) -> None
        """A DECIDE record: what was being decided, the evidence, the answer."""
        if isinstance(evidence, dict):
            evidence_text = "; ".join(
                "%s %s" % (name, _format_value(value)) for name, value in evidence.items()
            )
        elif isinstance(evidence, (list, tuple)):
            evidence_text = "; ".join(str(item) for item in evidence)
        else:
            evidence_text = str(evidence)
        evidence_text = evidence_text.strip()
        self._decisions.append(
            {"what": what, "evidence": evidence_text, "decision": decision}
        )
        if evidence_text:
            self.record("OBSERVED", what, evidence=evidence_text)
        self.record("DECIDE", what, value=decision)

    # -- QC actions --------------------------------------------------------
    def qc(self, check_name, plain_english, before, after, reason=None,
           changed=None, matched=None, outcome=None, warn=False, step=None):
        """The ONLY way to record a QC action.

        `before` and `after` are required positional arguments; that is what
        makes rule 5 of the plan structurally enforced rather than a convention.
        Calling qc() without them is a TypeError, not a silent omission.

        A zero-effect check becomes one PASS line.  A check that changed data
        records the observation, action, and resulting row counts separately.
        """
        if before is None or after is None:
            raise TypeError(
                "logger.qc() needs the variant count before and after the check; "
                "a QC action without counts cannot be told apart from one that never ran"
            )
        before = int(before)
        after = int(after)
        removed = before - after

        changed_count = int(changed or 0)
        matched_count = int(matched or 0)

        record = {
            "step": step or (self._current_step.label if self._current_step else None),
            "check": check_name,
            "reason": reason,
            "before": before,
            "after": after,
            "removed": removed,
            "changed": changed_count,
            "matched": matched_count,
            "outcome": outcome,
        }
        self._qc_actions.append(record)

        if removed == 0 and changed_count == 0 and matched_count == 0 and not warn:
            self.record("PASS", check_name, affected=0, rows=after)
            return record

        if removed < 0:
            self.record(
                "WARNING", check_name, issue="row_count_increased",
                rows_in=before, rows_out=after, added=-removed,
            )
        else:
            marker = "WARNING" if warn else "OBSERVED"
            self.record(
                marker, check_name,
                affected=max(removed, changed_count, matched_count),
                removed=max(0, removed), changed=changed_count,
                matched=matched_count,
            )
        action = str(plain_english).strip()
        if removed > 0:
            self.record(
                "ACTION", "remove_variants", count=removed,
                reason=reason or check_name,
            )
        elif matched_count > 0 and outcome == "kept":
            self.record(
                "ACTION", "keep_variants", count=matched_count,
                reason=reason or check_name,
            )
        elif matched_count > 0 and outcome == "fail":
            self.record(
                "ACTION", "fail_check", count=matched_count,
                reason=reason or check_name,
            )
        elif action:
            self.log("ACTION", action)
        self.record(
            "RESULT", check_name, rows_in=before, rows_out=after,
            removed=max(0, removed), changed=changed_count,
            matched=matched_count, outcome=outcome, reason=reason,
        )
        return record

    # -- steps -------------------------------------------------------------
    @contextmanager
    def step(self, number, total, title, func_name, rows_in=None, policy_keys=None):
        """Wrap one pipeline step.

        Logs STEP, INPUT and PARAM records; yields a small mutable context;
        then logs STATUS with in/out/removed/elapsed, or, if the
        step raised, FAILED with a plain-English explanation, the Python error
        and the traceback - and re-raises.
        """
        ctx = StepContext(self, number, total, title, func_name, rows_in)
        previous_step = self._current_step
        self._current_step = ctx
        label = _step_number(number, total)

        self.record("STEP", label, name=title, function=func_name)
        if rows_in is not None:
            self.record("INPUT", "variants", rows=rows_in)
        if policy_keys:
            self.settings(policy_keys)

        record = {
            "number": number,
            "total": total,
            "title": title,
            "func": func_name,
            "rows_in": rows_in,
            "rows_out": None,
            "removed": None,
            "elapsed": None,
            "status": "running",
        }
        self._steps.append(record)
        self._progress_event("start_step", number, total, title)

        try:
            yield ctx
        except BaseException as exc:
            elapsed = time.time() - ctx.started
            record["elapsed"] = elapsed
            record["status"] = "failed"
            record["error"] = "%s: %s" % (type(exc).__name__, exc)
            what, todo = _explain_exception(exc)
            if ctx.failure_hint:
                what = str(ctx.failure_hint)
            self.record("STATUS", label, status="FAILED", duration=_elapsed(elapsed))
            self.log("FAILED", "%s: %s" % (label, what))
            self.log("ACTION", todo, level="ERROR")
            self.log(
                "",
                "For developers - Python error: %s: %s" % (type(exc).__name__, exc),
                indent=1, level="ERROR", screen=False,
            )
            self.log("", "For developers - traceback:", indent=1, level="ERROR", screen=False)
            trace = traceback.format_exc().rstrip()
            for line in trace.split("\n"):
                self.log("", line, indent=2, level="ERROR", wrap=False, screen=False)
            self._progress_event("fail_step", number, total, title)
            self._current_step = previous_step
            raise
        else:
            elapsed = time.time() - ctx.started
            rows_out = ctx.rows_out if ctx.rows_out is not None else ctx.rows_in
            removed = ctx.removed
            if removed is None and ctx.rows_in is not None and rows_out is not None:
                removed = ctx.rows_in - rows_out
            record["rows_out"] = rows_out
            record["removed"] = removed
            record["elapsed"] = elapsed
            record["status"] = "ok"
            record["extra"] = dict(ctx.extra)

            outcome = ctx.extra.get("outcome")
            if outcome:
                self.record("RESULT", "stage_outcome", step=label, **(ctx.logged_outcome or outcome))

            self.record(
                "STATUS", label, status="OK", rows_in=ctx.rows_in,
                rows_out=rows_out, removed=removed, duration=_elapsed(elapsed),
            )
            self._progress_event(
                "complete_step",
                number,
                total,
                title,
                outcome["message"] if outcome else None,
                ctx.outcome_fields,
            )
            self._current_step = previous_step

    # -- output ------------------------------------------------------------
    def screen_text(self):
        # type: () -> str
        """Everything this worker would have printed, for the parent to print as
        one uninterrupted block."""
        return "\n".join(self._screen)

    def summary(self):
        # type: () -> Dict[str, Any]
        """Counts and timings, for the parent's footer and the run manifest."""
        steps = [dict(record) for record in self._steps]
        return {
            "sample_id": self.sample_id,
            "scope": self.scope,
            "label": self.label,
            "log_path": self.log_path,
            "warnings": self.warning_count,
            "errors": self.error_count,
            "steps": steps,
            "step_timings": dict(
                ("%s %s" % (_step_number(s["number"], s["total"]), s["title"]), s["elapsed"])
                for s in steps
            ),
            "qc_actions": [dict(record) for record in self._qc_actions],
            "decisions": [dict(record) for record in self._decisions],
            "elapsed": time.time() - self.started,
            "failed": any(s["status"] == "failed" for s in steps) or self.error_count > 0,
        }

    def close(self):
        """Safe to call twice, and safe to call after a failure."""
        self._progress_event("close")
        handle, self._fh = self._fh, None
        if handle is not None and not handle.closed:
            try:
                handle.flush()
                handle.close()
            except Exception:
                pass
        self._closed = True

    def __repr__(self):
        return "<PipelineLogger %s %s -> %s>" % (self.sample_id, self.label, self.log_path)


# =============================================================================
# Parent-side screen block
# =============================================================================


def _integer(value):
    # type: (Any) -> int
    try:
        return int(value or 0)
    except (TypeError, ValueError):
        return 0


def _sum_mapping(values):
    # type: (Any) -> int
    return sum(_integer(value) for value in (values or {}).values())


def _step_extra(summary, number):
    # type: (Dict[str, Any], int) -> Dict[str, Any]
    for step in summary.get("steps") or []:
        if _integer(step.get("number")) == number:
            return dict(step.get("extra") or {})
    return {}


def _source_name(value):
    # type: (Any) -> str
    text = str(value or "not available")
    if text in ("missing", "not available", "NA"):
        return "not provided"
    if text.startswith("column:"):
        return "column %s" % text.split(":", 1)[1]
    if text.startswith("fixed_value:"):
        return "fixed value %s" % text.split(":", 1)[1]
    if text.startswith("study_column:"):
        return "study column %s" % text.split(":", 1)[1]
    if text.startswith("external_file:"):
        path = text.split(":", 1)[1]
        return "external file %s" % os.path.basename(path) if path else "external file"
    if text == "calculated_from_beta_se":
        return "calculated as BETA / SE"
    if text == "built_from_coordinates":
        return "built from chromosome, position and alleles"
    return text.replace("_", " ")


def _analysis_line(label, detail, width, kind="analysis"):
    # type: (str, str, int, str) -> List[str]
    return screen_field(
        kind, label, detail, width=width, indent=6, label_width=20,
    ).splitlines()


def _chromosome_section(title, fields, width, kind="analysis"):
    """Render recorded stage metrics without parsing prose or scanning data."""
    lines = ["", screen_line(kind, title, indent=8)]
    for label, value in fields:
        lines.extend(screen_field(
            "info", label, str(value), width=width, indent=12, label_width=26,
        ).splitlines())
    return lines


def _chromosome_analysis_lines(summary, width):
    """Show recorded chromosome stages in execution order."""
    stages = summary.get("stage_qc") or {}
    lines = []

    def section(title, fields, kind="analysis"):
        lines.extend(_chromosome_section(title, fields, width, kind))

    def value(record, key):
        return record.get(key) if record.get(key) is not None else "not recorded"

    def count(record, key):
        return _count(record[key]) if record.get(key) is not None else "not recorded"

    effect = stages.get("beta_or_oddsratio_qc") or {}
    if effect:
        is_or = effect.get("effect_type") == "odds_ratio"
        section("Effect-scale conversion", [
            ("Input column", value(effect, "effect_col")),
            ("Input effect type", value(effect, "effect_type")),
            ("Decision source", value(effect, "effect_decision_source")),
            ("Transformation", "No input effect supplied; derivation occurs later" if effect.get("status") == "skipped_no_input_effect_column"
             else "BETA = log(OR)" if is_or else "BETA used directly"),
            ("Rows processed", count(effect, "initial_variants")),
            ("Standard error", "Converted from OR scale to log-odds" if effect.get("se_rescaled_by_or")
             else "%s scale; unchanged" % value(effect, "se_input_scale")),
            ("Non-positive OR detected", count(effect, "or_non_positive_count") if is_or else "Not applicable"),
            ("Non-positive OR action", value(effect, "or_non_positive_action") if is_or else "Not applicable"),
            ("Retained at this step", count(effect, "final_total")),
        ])

    eaf = stages.get("eaf_qc") or {}
    strand = eaf.get("strand_orientation") or {}
    if strand:
        actions = strand.get("actions") or {}
        section("Strand and allele alignment", [
            ("Reference", os.path.basename(str(strand.get("reference_file") or "not recorded"))),
            ("Reference population", value(strand, "reference_population_column")),
            ("Mode / consensus", "%s / %s" % (value(strand, "strand_mode"), value(strand, "study_strand_consensus"))),
            ("Matched", _count(sum(_integer(v) for v in actions.values()))
             + (" / %s" % _count(strand["initial_variants"]) if strand.get("initial_variants") is not None else "")),
            ("Alleles unchanged", _count(_integer(actions.get("forward")))),
            ("Alleles swapped", _count(_integer(actions.get("forward_swapped")))),
            ("Reverse-complement", _count(_integer(actions.get("reverse_complement")))),
            ("Complement + swap", _count(_integer(actions.get("reverse_complement_swapped")))),
            ("Unmatched removed", count(strand, "reference_unmatched")),
            ("Unmatched retained", count(strand, "reference_unmatched_retained")),
            ("Palindromic conflicts", "%s removed; study-wide strand and allele-frequency evidence supported opposite orientations" % _count(_integer(strand.get("palindromic_frequency_conflict")))),
            ("Orientation unavailable", _count(_integer(strand.get("palindromic_orientation_unavailable")))),
            ("Palindromic AF failures", _count(_integer(strand.get("palindromic_frequency_discordant")))),
            ("Palindromic ambiguous", count(strand, "palindromic_ambiguous")),
            ("Reference ambiguous", count(strand, "reference_ambiguous")),
            ("Retained", count(strand, "final_variants")),
        ], "genetic")

    if eaf:
        external = str(eaf.get("decision_source") or "").startswith("external")
        interpretation = "EAF" if external or eaf.get("maf_reference_decision") == "eaf" or eaf.get("study_decision_eaf_is_maf") is False else "MAF suspected / unresolved"
        fields = [
            ("Source", eaf.get("external_eaf_file") if external else "Study column %s" % value(eaf, "final_eaf_col")),
            ("Frequency interpretation", interpretation),
        ]
        if external:
            fields += [
                ("Reference column", value(eaf, "external_eaf_column")),
                ("Direct matches", count(eaf, "external_merge_direct_match_rows")),
                ("Swapped matches", count(eaf, "external_merge_flip_match_rows")),
                ("Unmatched", count(eaf, "external_merge_unmatched_rows")),
            ]
        else:
            actions = strand.get("actions") or {}
            swapped = _integer(actions.get("forward_swapped")) + _integer(actions.get("reverse_complement_swapped"))
            fields.append(("Transformation", "%s allele-swapped frequencies inverted using 1 − EAF" % _count(swapped)
                           if interpretation == "EAF" else "No EAF inversion asserted while frequency type is unresolved"))
            outside = eaf.get("internal_af_outside_unit_interval_count", eaf.get("internal_af_out_of_range_count"))
            fields.append(("Out-of-range detected", _count(outside) if outside is not None else "not recorded"))
        if eaf.get("strand_af_tolerance") is not None:
            fields.append(("Reference AF tolerance", "|aligned study EAF − reference AF| ≤ %s" % eaf["strand_af_tolerance"]))
        for label, prefix, action_key in [
            ("Non-palindromic AF", "strand_non_palindromic_af", "strand_af_discordance_action"),
            ("Palindromic AF", "strand_palindromic_af", "strand_palindromic_af_discordance_action"),
        ]:
            if eaf.get(prefix + "_discordant") is not None:
                action = eaf.get(action_key)
                handling = {"warn": "retained with warning; no frequencies changed by this check",
                            "reject": "removed", "fail": "stop on disagreement"}.get(action, "action %s" % action)
                fields.append((label, "%s / %s comparable variants disagree; %s" % (
                    count(eaf, prefix + "_discordant"), count(eaf, prefix + "_comparable"), handling)))
        fields += [
            ("Total strand/EAF removals", count(eaf, "variants_removed")),
            ("Retained", count(eaf, "final_variants")),
        ]
        if strand.get("final_variants") is not None and eaf.get("final_variants") is not None:
            fields.insert(-1, ("Additional AF removals", _count(strand["final_variants"] - eaf["final_variants"])))
        section("Allele frequency", fields, "genetic")

    size = stages.get("sample_size_qc") or {}
    if size:
        neff = str(size.get("Neff_status") or "not recorded")
        if neff == "calculated_from_case_control":
            neff = "NEFF = 4 / (1/Ncase + 1/Ncontrol)"
        elif neff.startswith("fallback_ncontrol_only"):
            neff = "Control-count input used as total sample size"
        elif neff.startswith("fallback_ncase_only"):
            neff = "Case-count input used as total sample size"
        section("Sample size", [
            ("Cases", _source_name(size.get("ncase_source"))),
            ("Controls / total N", _source_name(size.get("ncontrol_source"))),
            ("Calculation", neff), ("Minimum count", value(size, "min_value")),
            ("Removed", _count(_sum_mapping(size.get("removed_by_reason")))),
        ], "count")

    inputs = stages.get("effect_from_z_qc") or {}
    if inputs:
        if inputs.get("beta_computed") and inputs.get("se_computed"):
            handling = "standardized BETA estimate and SE derived from Z, EAF and NEFF"
            if inputs.get("effect_estimate_citation"):
                handling += "; " + inputs["effect_estimate_citation"]
        elif inputs.get("beta_computed"):
            handling = "BETA derived as Z × SE"
        elif inputs.get("se_computed") or inputs.get("se_recovered_from_z"):
            handling = "positive SE derived as abs(BETA / Z) after signed agreement"
        elif inputs.get("has_beta") and inputs.get("has_se"):
            handling = "Using available BETA and SE; nothing derived from Z"
        else:
            handling = "SE deferred to the p-value stage"
        section("Effect inputs", [
            ("Effect estimate", "%s%s" % (value(inputs, "effect_column"),
             " (BETA already converted from odds ratios)" if effect.get("effect_type") == "odds_ratio" else "")),
            ("Standard error column", value(inputs, "se_column")),
            ("Handling", handling),
            ("SE recovered from Z", count(inputs, "se_recovered_from_z")),
            ("BETA/Z sign mismatches", count(inputs, "beta_z_sign_mismatches")),
            ("Sign mismatch action", value(inputs, "beta_z_sign_mismatch_action")),
            ("Removed", count(inputs, "variants_removed_total")),
        ])

    pvalue = stages.get("pval_detection_and_conversion_qc") or {}
    if pvalue:
        section("P-values", [
            ("Input column / scale", "%s / %s" % (value(pvalue, "pvalue_column"), value(pvalue, "detected_scale"))),
            ("Decision source", value(pvalue, "pvalue_decision_source")),
            ("Accepted range", "%s to %s" % (value(pvalue, "pvalue_clip_low"), value(pvalue, "pvalue_clip_high"))),
            ("Outside-range action", value(pvalue, "pvalue_out_of_range_action")),
            ("Clipped", _count(_integer(pvalue.get("variants_with_pvalues_clipped_low")) + _integer(pvalue.get("variants_with_pvalues_clipped_high")))),
            ("Removed", _count(max(0, _integer(pvalue.get("initial_total_variants_with_pvalues")) - _integer(pvalue.get("after_filter_variants_with_pvalues"))))),
        ])

    se = stages.get("se_from_beta_pvalue_qc") or {}
    if se:
        status = se.get("status")
        handling = "supplied SE used directly" if status == "SE already present" else (
            "Calculated from BETA and %s-sided p-values" % value(se, "se_tail") if status == "SE calculated" else status or "not recorded")
        section("Standard error", [
            ("Handling", handling), ("Undefined", count(se, "se_undefined")),
            ("Removed", _count(max(0, _integer(se.get("initial_variants")) - _integer(se.get("total_variants"))))),
        ])

    validation = stages.get("effect_statistics_validation_qc") or {}
    zscore = stages.get("z_from_beta_se_qc") or {}
    if zscore:
        section("Z score", [
            ("Calculation / source", _source_name(zscore.get("z_source"))),
            ("Exact zero effects", count(validation, "variants_with_zero_beta")),
            ("Zero-effect action", value(validation, "beta_zero_action")),
            ("Invalid BETA/SE removed", count(validation, "variants_removed_invalid_beta_se")),
        ])
    if validation:
        concordance = validation.get("z_pval_concordance") or {}
        section("Effect-statistic QC", [
            ("Checked", "BETA, SE and Z"), ("Z/P disagreement action", value(concordance, "action")),
            ("Z/P comparable", count(concordance, "checked")),
            ("Z/P disagreements", count(concordance, "discordant")),
            ("Removed", count(validation, "removed_total")),
        ])

    info = stages.get("info_qc") or {}
    if info:
        mach = info.get("info_score_type") == "mach_rsq"
        section("Imputation quality", [
            ("Source", _source_name(info.get("source"))),
            ("Input column", value(info, "info_column")),
            ("Score type", value(info, "info_score_type")),
            ("Accepted range", "%s–%s" % (value(info, "info_clip_min"), value(info, "info_mach_rsq_max" if mach else "info_clip_max"))),
            ("Tolerance corrections", "Not applied to MaCH Rsq" if mach else "%s values above %s and ≤%s set to %s" % (
                count(info, "rescaled_within_tolerance"), value(info, "info_clip_max"),
                value(info, "info_clip_tolerance"), value(info, "info_clip_max"))),
            ("Outside-range action", value(info, "info_out_of_range_action")),
            ("Outside-range detected", count(info, "out_of_range")),
            ("Missing / action", "%s / %s" % (count(info, "missing_info_after"), value(info, "info_on_missing_action"))),
            ("Removed", _count(_integer(info.get("rejected_out_of_range")) + _integer(info.get("rejected_missing")))),
        ])

    identifiers = _step_extra(summary, 12)
    if identifiers:
        filled = identifiers.get("identifiers_filled")
        if filled is None and identifiers.get("missing_identifiers") == 0:
            filled = 0
        section("Variant IDs", [
            ("Source", _source_name(identifiers.get("source"))),
            ("Separators rewritten", count(identifiers, "identifiers_rewritten")),
            ("Missing IDs filled", _count(filled) if filled is not None else "not recorded"),
        ], "genetic")
    final = stages.get("final_completeness_qc") or {}
    if final:
        section("Final completeness", [
            ("Required fields", ", ".join(final.get("required_fields") or [])),
            ("Missing observations", _count(_sum_mapping(final.get("missing_by_field")))),
            ("Missing-value action", value(final, "on_missing")),
            ("Removed", count(final, "removed_total")),
        ])
    return lines


def _target_vcf(summary):
    # type: (Dict[str, Any]) -> Optional[str]
    for step in summary.get("steps") or []:
        candidate = (step.get("extra") or {}).get("target_vcf")
        if candidate:
            return str(candidate)
    return None


def _liftover_counts(summary):
    # type: (Dict[str, Any]) -> Tuple[Optional[int], Optional[int]]
    for action in summary.get("qc_actions") or []:
        if str(action.get("check") or "").strip().upper() == "LIFTED":
            return int(action.get("before") or 0), int(action.get("after") or 0)
    return None, None


def _screen_messages(screen_text):
    # type: (str) -> List[str]
    messages = []
    for line in (screen_text or "").splitlines():
        cleaned = line[_PREFIX_WIDTH:].strip() if line.startswith("[") else line.strip()
        if cleaned and cleaned not in messages:
            messages.append(cleaned)
    return messages


def print_chromosome_summary(
    chromosome,
    sample_id,
    attempt,
    status,
    elapsed,
    screen_text,
    summary,
    stream=None,
    width=None,
):
    # type: (Any, str, Any, str, Optional[float], str, Optional[Dict[str, Any]], Any, Optional[int]) -> str
    """Print one concise, user-facing chromosome outcome.

    Only the parent calls this, from its single-threaded loop, which is what
    makes interleaving between chromosomes impossible. Detailed steps remain
    in the chromosome log. Returns the text for callers that buffer blocks.
    """
    if stream is None:
        stream = sys.stdout
    screen_width = int(width or _LINE_WIDTH)

    status_text = str(status).upper()
    chromosome_name = _scope_label(chromosome).replace("chr", "Chromosome ", 1)
    lines = [
        screen_line("genetic", chromosome_name, indent=4),
        screen_field(
            "success" if status_text == "OK" else "error",
            "Status",
            "%s in %s (attempt %s)"
            % ("completed" if status_text == "OK" else "failed", _elapsed(elapsed), attempt),
            width=screen_width,
            indent=6,
            label_width=20,
        ),
    ]
    summary = summary or {}
    rows_in = summary.get("rows_in")
    rows_out = summary.get("rows_out")
    rejected = summary.get("rejected")
    if rows_in is not None:
        lines.extend(_analysis_line(
            "Variants read", _count(rows_in), screen_width, "count",
        ))

    stage_lines = _chromosome_analysis_lines(summary, screen_width)
    if stage_lines:
        lines.append("")
        lines.extend(stage_lines)

    if rows_out is not None and rejected is not None:
        fields = [("Input", _count(rows_in)), ("Total removed", _count(rejected))]
        fields.extend((reason.replace("_", " "), _count(count))
                      for reason, count in (summary.get("reject_counts") or {}).items() if count)
        fields.append(("Harmonised", _count(rows_out)))
        lines.extend(_chromosome_section("Harmonisation accounting", fields, screen_width, "count"))

    lifted_in, lifted_out = _liftover_counts(summary)
    accounting = (_step_extra(summary, 16) or {}).get("liftover_accounting") or {}
    if lifted_in is not None or accounting:
        fields = [("Sent for conversion", _count(rows_out)),
                  ("Target VCF retained", _count(accounting.get("final", lifted_out)))]
        if accounting:
            fields.extend([
                ("Liftover input", _count(accounting.get("input"))),
                ("Liftover rejected", _count(accounting.get("rejected"))),
                ("Excluded after liftover", "%s successfully lifted variants with swapped alleles" % _count(accounting.get("swap_excluded"))),
                ("Exclusion policy", "vcf.liftover_swap = %s" % accounting.get("swap_policy")),
            ])
        elif lifted_in is not None:
            fields.append(("Net reduction", "%s; cause breakdown not recorded" % _count(lifted_in - lifted_out)))
        lines.extend(_chromosome_section("VCF creation and liftover", fields, screen_width, "genetic"))
    elif _target_vcf(summary):
        lines.extend(_chromosome_section("VCF creation and liftover", [
            ("Result", "Created and annotated; count breakdown not recorded"),
        ], screen_width, "genetic"))

    warnings = int(summary.get("warnings") or 0)
    errors = int(summary.get("errors") or 0)
    failed = status_text != "OK"
    if warnings or errors or failed:
        if failed and errors == 0:
            errors = 1
        lines.extend(_analysis_line("Messages", "%s warning%s, %s error%s; see the log." % (
            _count(warnings), "" if warnings == 1 else "s",
            _count(errors), "" if errors == 1 else "s",
        ), screen_width, "error" if errors else "warning"))
        af_counter_continuation = False
        for message in _screen_messages(screen_text):
            # These retained-AF warnings have a complete labelled breakdown
            # above. Do not repeat the generic logger counters as prose.
            if not failed and (summary.get("stage_qc") or {}).get("eaf_qc"):
                if "study/reference frequency concordance" in message:
                    af_counter_continuation = True
                    continue
                if af_counter_continuation and re.fullmatch(r"(?:affected|removed|changed|matched)=\d+(?:\s+(?:affected|removed|changed|matched)=\d+)*", message):
                    continue
            af_counter_continuation = False
            lines.extend(_analysis_line(
                "Reason", message, screen_width,
                "error" if errors else "warning",
            ))
    else:
        lines.extend(_analysis_line(
            "Messages", "no warnings or errors.", screen_width, "info",
        ))

    target = _target_vcf(summary)
    file_names = []
    if target:
        file_names.append(("Target VCF", os.path.basename(target)))
    log_path = summary.get("log_path")
    file_names.append(("Detailed log", os.path.basename(log_path) if log_path else "not recorded"))
    lines.extend(_chromosome_section("Files", file_names, screen_width, "info"))

    text = "\n".join(lines)
    stream.write("\n" + text + "\n\n")
    try:
        stream.flush()
    except Exception:
        pass
    return text
