"""Shared execution plumbing for harmonisation steps.

This module deliberately contains no scientific thresholds or formulas.
"""

from __future__ import annotations

from contextlib import contextmanager
from typing import Any, Iterable

import polars as pl

from postgwas.core.ui.screen import screen_field, screen_line
from postgwas.core.ui.progress import print_screen_block
from postgwas.core.values import format_count, optional_text

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


def active_context(logger=None, existing=None, rows_in=None):
    """Return the active step, its logger, or a silent compatible context.

    Scientific helpers are also called directly in focused tests and library
    use, outside a ``logger.step(...)`` block.  This resolver gives those calls
    the same message/QC interface without duplicating a module-specific sink.
    It does not open or close a step; :func:`step_context` owns that lifecycle.
    """
    if existing is not None:
        return existing
    if logger is not None:
        return logger
    return NullStepContext(rows_in)


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
            print_screen_block(screen_line(
                "warning" if warn else "info",
                message,
                indent=4 + indent * 2,
            ))
        return
    method = logger.warn if warn else logger.info
    method(message, indent=indent)


def emit_high_visibility_warning(logger, message, *, screen_message=None) -> None:
    """Record a plain warning and display it through the shared warning theme."""
    if logger is not None:
        try:
            logger.log("WARNING", message, screen=False)
        except (AttributeError, TypeError):
            logger.warn(message)
    rendered = screen_message if screen_message is not None else screen_line("warning", message, indent=4)
    print_screen_block(rendered)


def format_inference_warning(
    title, *, input_column, inferred_type, evidence, planned_action, reason,
    recommendation,
) -> str:
    """Render study-inference evidence using one shared screen hierarchy.

    Callers supply already-recorded evidence and scientific explanations. This
    formatter applies no thresholds, makes no decisions and reads no data.
    """
    lines = [screen_line("warning", title, indent=4)]
    for label, value in (("Input column", input_column), ("Inferred type", inferred_type)):
        lines.append(screen_field("info", label, str(value), indent=8, label_width=20))
    lines.append(screen_line("analysis", "Evidence used", indent=8))
    for label, value in evidence:
        lines.append(screen_field("info", label, str(value), indent=12, label_width=20))
    for label, value in (
        ("Planned action", planned_action), ("Why this warning", reason),
        ("Recommended", recommendation),
        ("Run behaviour", "Continuing with the inferred type."),
    ):
        lines.append(screen_field("info", label, value, indent=8, label_width=20))
    return "\n".join(lines)


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
    condition: pl.Expr | pl.Series,
    *,
    step_label: str,
    reason: str,
    context=None,
    collector=None,
    counters: dict[str, int] | None = None,
    detail: str | None = None,
    check_name: str | None = None,
    description: str | None = None,
    warn_on_remove: bool = False,
    warn_without_collector: bool = False,
    record_empty: bool = False,
) -> tuple[pl.DataFrame, int]:
    """Reject matching rows once, log the action, and optionally update counters.

    A null rejection predicate is treated as a match, preserving the fail-closed
    contract of :class:`RejectCollector`.  Domain rules where null means "this
    rule does not apply" must make that scientific decision explicitly with
    ``condition.fill_null(False)`` before calling this function.  This keeps a
    missing frequency, for example, out of the out-of-range category so the
    later missing-frequency policy can handle it with the correct provenance.
    """
    before = frame.height
    if before == 0 and not record_empty:
        return frame, 0
    if collector is not None:
        kept = collector.reject(
            frame, condition, step_label, reason, detail=detail
        )
    else:
        kept = frame.filter(~condition.fill_null(True))
    removed = before - kept.height
    collector_logged = (
        collector is not None and getattr(collector, "logger", None) is not None
    )
    if context is not None and not collector_logged:
        explanation = description or REASONS[reason]
        if detail and description is None:
            explanation = "%s (%s)" % (explanation, detail)
        context.qc(
            check_name or reason.replace("_", " "),
            explanation,
            before,
            kept.height,
            reason=reason,
            warn=bool(warn_on_remove and removed > 0),
            step=step_label,
        )
    if (
        collector is None
        and warn_without_collector
        and removed > 0
        and context is not None
    ):
        context.warn(
            "No rejected-variants collector was supplied to this step, so these "
            "%s variants were removed without being written to the "
            "rejected-variants file." % format_count(removed)
        )
    if counters is not None:
        counters[reason] = counters.get(reason, 0) + removed
    return kept, removed


__all__ = [
    "NullLogger",
    "NullStepContext",
    "active_context",
    "configured_column",
    "emit_message",
    "emit_high_visibility_warning",
    "format_inference_warning",
    "log_info",
    "log_warning",
    "reject_rows",
    "resolve_policies",
    "step_context",
]
