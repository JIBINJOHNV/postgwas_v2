"""Shared contracts for pipeline validation performed before execution."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any

from postgwas.core.input_validation import (
    PreflightFileIdentity,
    capture_preflight_file_identities,
    require_unchanged_preflight_files,
)
from postgwas.core.vcf import IndexedVcfValidation


@dataclass(frozen=True)
class PipelinePreflightEvidence:
    """Validated entry data, external resources, and deferred handoff checks."""

    module: str
    input_vcf: Mapping[str, Any]
    resources: Any = field(default_factory=dict)
    deferred_checks: tuple[str, ...] = ()


@dataclass(frozen=True)
class PreflightLogEvent:
    """One structured log call retained until the canonical log is available."""

    method: str
    arguments: tuple[Any, ...]
    fields: tuple[tuple[str, Any], ...] = ()


class PreflightLogBuffer:
    """Capture validator logging without creating analysis output directories."""

    def __init__(self) -> None:
        self._events: list[PreflightLogEvent] = []

    @property
    def events(self) -> tuple[PreflightLogEvent, ...]:
        return tuple(self._events)

    def record(self, marker, subject, **values) -> None:
        self._events.append(PreflightLogEvent(
            "record", (marker, subject), tuple(values.items()),
        ))

    def info(self, message, indent=0) -> None:
        self._events.append(PreflightLogEvent(
            "info", (message,), (("indent", indent),),
        ))

    def warn(self, message, indent=0) -> None:
        self._events.append(PreflightLogEvent(
            "warn", (message,), (("indent", indent),),
        ))

    warning = warn


def replay_preflight_log(events, logger) -> None:
    """Replay buffered validation events into the module's canonical log."""
    for event in events:
        method = getattr(logger, event.method)
        method(*event.arguments, **dict(event.fields))


def require_pipeline_input_vcf(
    preflight_evidence: Mapping[str, Any] | None,
) -> Mapping[str, Any]:
    """Return the common entry-VCF evidence or reject an invalid call order."""
    evidence = preflight_evidence or {}
    validated = evidence.get("input_vcf")
    if not isinstance(validated, Mapping):
        raise RuntimeError(
            "Pipeline module preflight requires the shared input-VCF "
            "validation to complete first."
        )
    indexed = validated.get("indexed")
    harmonised = validated.get("harmonised")
    if not isinstance(indexed, IndexedVcfValidation) or not isinstance(
        harmonised, Mapping
    ):
        raise RuntimeError(
            "Shared input-VCF evidence is incomplete; indexed and harmonised "
            "contracts must both be validated before module resources."
        )
    return validated


def pipeline_preflight_evidence(
    module: str,
    preflight_evidence: Mapping[str, Any] | None,
    *,
    resources: Any = None,
    deferred_checks: tuple[str, ...] = (),
) -> PipelinePreflightEvidence:
    """Build the uniform result returned by every module preflight callable."""
    name = str(module).strip()
    if not name:
        raise ValueError("Pipeline preflight evidence requires a module name")
    return PipelinePreflightEvidence(
        module=name,
        input_vcf=require_pipeline_input_vcf(preflight_evidence),
        resources={} if resources is None else resources,
        deferred_checks=tuple(str(check) for check in deferred_checks),
    )


__all__ = [
    "PipelinePreflightEvidence",
    "PreflightFileIdentity",
    "PreflightLogBuffer",
    "PreflightLogEvent",
    "capture_preflight_file_identities",
    "pipeline_preflight_evidence",
    "replay_preflight_log",
    "require_unchanged_preflight_files",
    "require_pipeline_input_vcf",
]
