"""Shared execution plumbing for harmonisation steps.

This module deliberately contains no scientific thresholds or formulas.
"""

from __future__ import annotations

import sys
from contextlib import contextmanager
from typing import Any, Iterable

import polars as pl

from postgwas.core.ui.screen import screen_line
from postgwas.core.values import optional_text

from ..policies import default_policies
from ..rejects import REASONS


class NullStepContext:
    """StepContext-compatible sink used by direct, logger-free function calls."""

    def __init__(self, rows_in=None):
        self.rows_in = rows_in
        self.rows_out = rows_in
        self.removed = 0
        self.extra = {}
        self.failure_hint = None

    def info(self, message, indent=0):
        return None

    debug = info

    def warn(self, message, indent=0):
        return None

    warning = warn
    error = warn
    skip = warn

    def decide(self, what, evidence, decision):
        return None

    def input(self, subject, **values):
        return None

    observed = input
    action = input
    output = input

    def qc(self, check_name, plain_english, before, after, **kwargs):
        if before is None or after is None:
            raise TypeError(
                "QC check %r requires variant counts before and after." % check_name
            )
        return {
            "check": check_name,
            "before": int(before),
            "after": int(after),
            "removed": int(before) - int(after),
        }

    def set_rows(self, rows_out, removed=None):
        self.rows_out = rows_out
        self.removed = (
            self.rows_in - rows_out
            if removed is None and self.rows_in is not None and rows_out is not None
            else removed
        )


class NullLogger:
    """PipelineLogger-compatible sink for direct scientific function calls."""

    @contextmanager
    def step(self, number, total, title, operation, rows_in=None, policy_keys=None):
        yield NullStepContext(rows_in)

    def info(self, message, indent=0):
        return None

    debug = info

    def warn(self, message, indent=0):
        return None

    warning = warn
    error = warn
    skip = warn

    def decide(self, what, evidence, decision):
        return None

    def record(self, marker, subject, **values):
        return None

    def qc(self, check_name, plain_english, before, after, **kwargs):
        return NullStepContext().qc(
            check_name, plain_english, before, after, **kwargs
        )


def resolve_policies(policies=None):
    """Use the supplied immutable policies or load canonical defaults."""
    return policies if policies is not None else default_policies()


def log_info(logger, message) -> None:
    """Record an informational message when a logger is present."""
    if logger is not None:
        logger.info(message)


def log_warning(logger, message) -> None:
    """Record a warning when a logger is present."""
    if logger is not None:
        logger.warn(message)


def emit_message(
    logger,
    message,
    warn: bool = False,
    indent: int = 0,
    *,
    screen_when_unlogged: bool = False,
) -> None:
    """Log one message, optionally showing it for a logger-free direct call."""
    if logger is None:
        if screen_when_unlogged:
            print(screen_line(
                "warning" if warn else "info",
                message,
                indent=4 + indent * 2,
            ))
        return
    method = logger.warn if warn else logger.info
    method(message, indent=indent)


def emit_high_visibility_warning(logger, message) -> None:
    """Record a plain warning and show it in bright red on an interactive terminal."""
    if logger is not None:
        try:
            logger.log("WARNING", message, screen=False)
        except (AttributeError, TypeError):
            logger.warn(message)
    rendered = screen_line("warning", message, indent=4)
    if getattr(sys.stdout, "isatty", lambda: False)():
        rendered = "\033[1;91m%s\033[0m" % rendered
    sys.stdout.write(rendered + "\n")
    sys.stdout.flush()


def configured_column(mapping: dict, keys: Iterable[str], frame: pl.DataFrame) -> str | None:
    """Return the first configured candidate that exists in the frame."""
    if not isinstance(mapping, dict):
        return None
    for key in keys:
        name = optional_text(mapping.get(key))
        if name and name in frame.columns:
            return name
    return None


@contextmanager
def step_context(
    logger,
    existing,
    *,
    number: int,
    total: int,
    title: str,
    operation: str,
    rows_in: int | None,
    policy_keys=None,
):
    """Use an existing step, open a logged step, or yield a silent step."""
    if existing is not None:
        yield existing
    elif logger is not None:
        with logger.step(
            number, total, title, operation,
            rows_in=rows_in, policy_keys=policy_keys,
        ) as created:
            yield created
    else:
        yield NullStepContext(rows_in)


def reject_rows(
    frame: pl.DataFrame,
    condition: pl.Expr,
    *,
    step_label: str,
    reason: str,
    context=None,
    collector=None,
    counters: dict[str, int] | None = None,
    detail: str | None = None,
    check_name: str | None = None,
    description: str | None = None,
) -> tuple[pl.DataFrame, int]:
    """Reject matching rows once, log the action, and optionally update counters."""
    before = frame.height
    if before == 0:
        return frame, 0
    if collector is not None:
        kept = collector.reject(
            frame, condition, step_label, reason, detail=detail
        )
    else:
        kept = frame.filter(~condition.fill_null(True))
        if context is not None:
            explanation = description or REASONS[reason]
            if detail and description is None:
                explanation = "%s (%s)" % (explanation, detail)
            context.qc(
                check_name or reason.replace("_", " "),
                explanation,
                before,
                kept.height,
                reason=reason,
                step=step_label,
            )
    removed = before - kept.height
    if counters is not None:
        counters[reason] = counters.get(reason, 0) + removed
    return kept, removed


__all__ = [
    "NullLogger",
    "NullStepContext",
    "configured_column",
    "emit_message",
    "emit_high_visibility_warning",
    "log_info",
    "log_warning",
    "reject_rows",
    "resolve_policies",
    "step_context",
]
