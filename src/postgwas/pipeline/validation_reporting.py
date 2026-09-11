"""One pipeline input-validation presentation and atomic startup audit."""

from __future__ import annotations

from pathlib import Path

from postgwas.core.preflight import PipelinePreflightEvidence
from postgwas.core.validation_reporting import validation_audit_path, write_validation_audit

# Audit ownership marker is a serialization protocol, not a user filename/policy.
_REPORT_KIND = "postgwas.pipeline.input_validation"


def validation_report_path(args, configuration) -> Path:
    """Keep the audit within the output root and away from the entry input."""
    path = validation_audit_path(args.output_directory, configuration.pipeline.validation.report_file)
    if getattr(args, "vcf", None) and path == Path(args.vcf).expanduser().resolve():
        raise ValueError("Pipeline validation report would overwrite the input VCF: %s" % path)
    return path


def save_pipeline_validation_report(
    session, evidence, modules, path, *, error=None, phase="pipeline_startup",
):
    """Persist only established startup evidence, including failure and deferral.

    A deferred item describes a later validation boundary. Its status is never
    promoted to passed merely because the pipeline or an unrelated check ran.
    """
    document = {"report_kind": _REPORT_KIND}
    stages = []
    for name in dict.fromkeys(modules):
        item = evidence.get(name)
        stages.append({
            "module": name,
            "preflight": "passed" if isinstance(item, PipelinePreflightEvidence) else "not_passed",
            "deferred_checks": list(item.deferred_checks) if isinstance(item, PipelinePreflightEvidence) else [],
        })
    document.update({
        "phase": phase,
        "status": "failed" if error is not None else "passed",
        "error": None if error is None else "%s: %s" % (type(error).__name__, error),
        "stages": stages,
        "interpretation": (
            "Each file lists only the checks performed. Availability or header/index "
            "checks do not establish all-record scientific validity. Declared build "
            "and population are not independently inferred from filenames. Generated "
            "inputs are validated at their consuming stages. The stages' deferred_checks "
            "describe startup deferrals, not a claim that those checks passed."
        ),
    })
    return write_validation_audit(document, session.records, path)


__all__ = ["validation_report_path", "save_pipeline_validation_report"]
