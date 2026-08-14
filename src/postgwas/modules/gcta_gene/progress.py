"""Measured progress from append-only native GCTA gene-test outputs."""

from __future__ import annotations

import re
from pathlib import Path

from postgwas.core.ui import MeasuredProgress


# mBAT-combo declares its mapped-gene total before testing and flushes one
# result row after each successfully analysed gene. Genes for which GCTA does
# not produce a valid result can be omitted, so the row count is a conservative
# lower bound while work is active. fastBAT instead logs its measured compute
# counter before writing the result table. These are upstream protocol
# invariants, not estimates made by PostGWAS.
_MBAT_TOTAL_PATTERNS = (
    re.compile(r"Running\s+mBAT-combo analysis for\s+(\d+)\s+gene\(s\)", re.I),
    re.compile(r"(\d+)\s+genes have been mapped to SNP data", re.I),
)
_FASTBAT_PROGRESS_PATTERN = re.compile(
    r"(\d+)\s+of\s+(\d+)\s+(?:genes|sets)\.", re.I,
)
_VARIANT_MATCH_PATTERN = re.compile(
    r"Matching the GWAS meta-analysis results to the genotype data", re.I,
)
_GENE_MAPPING_PATTERN = re.compile(
    r"Mapping the physical positions of genes to SNP data", re.I,
)
_FASTBAT_START_PATTERN = re.compile(r"Running\s+fastBAT analysis", re.I)
_MBAT_START_PATTERN = re.compile(r"Running\s+mBAT-combo analysis", re.I)


class _AppendedTextReader:
    """Read only bytes appended since the preceding sample."""

    def __init__(self, path: Path):
        self.path = Path(path)
        self.offset = 0
        self.text = ""

    def read(self) -> str:
        if not self.path.is_file():
            return ""
        size = self.path.stat().st_size
        if size < self.offset:
            self.offset = 0
            self.text = ""
        with self.path.open("rb") as handle:
            handle.seek(self.offset)
            appended = handle.read()
            self.offset = handle.tell()
        if appended:
            self.text += appended.decode("utf-8", errors="replace")
        return self.text


class _AppendedResultCounter:
    """Count complete nonempty rows without rescanning a growing result file."""

    def __init__(self, path: Path):
        self.path = Path(path)
        self.offset = 0
        self.partial = b""
        self.header_seen = False
        self.rows = 0

    def count(self) -> int:
        if not self.path.is_file():
            return self.rows
        size = self.path.stat().st_size
        if size < self.offset:
            self.offset = 0
            self.partial = b""
            self.header_seen = False
            self.rows = 0
        with self.path.open("rb") as handle:
            handle.seek(self.offset)
            appended = handle.read()
            self.offset = handle.tell()
        if not appended:
            return self.rows
        lines = (self.partial + appended).split(b"\n")
        self.partial = lines.pop()
        for line in lines:
            if not line.strip():
                continue
            if not self.header_seen:
                self.header_seen = True
            else:
                self.rows += 1
        return self.rows


class GctaResultProgress:
    """Translate GCTA's native plan and result rows into measured progress."""

    def __init__(
        self,
        *,
        method: str,
        unit_label: str,
        native_log: Path,
        result_file: Path,
        logger,
        enabled: bool,
        total_hint: int | None = None,
    ):
        self.method = str(method)
        self.unit_label = str(unit_label)
        self.logger = logger
        self._uses_result_rows = self.method == "mbat_combo"
        self.total = None if total_hint is None else int(total_hint)
        self._total_source = (
            (
                "gcta_native_mapped_plan"
                if self._uses_result_rows else "gcta_native_progress"
            )
            if total_hint is None else "validated_postgwas_set_plan"
        )
        if self.total is not None and self.total < 1:
            raise ValueError("GCTA progress total_hint must be a positive integer")
        self.completed = 0
        self._total_recorded = False
        self._native_total_mismatch_recorded = False
        self._last_logged_percentage = None
        self._native_log = _AppendedTextReader(native_log)
        self._results = (
            _AppendedResultCounter(result_file)
            if self._uses_result_rows else None
        )
        self._progress = MeasuredProgress(
            "GCTA %s measured progress" % self.method.replace("_", " "),
            enabled=enabled,
        )
        # These transitions are emitted by GCTA itself. They let PostGWAS name
        # indeterminate setup work without converting unequal phases into a
        # fabricated percentage. The native gene/set counter remains the only
        # source of measured completion.
        analysis_titles = {
            "fastbat_gene": "Evaluate annotated genes with fastBAT",
            "fastbat_segment": "Test genomic segments with fastBAT",
            "fastbat_set": "Test variant sets with fastBAT",
            "mbat_combo": "Test mapped genes",
        }
        analysis_pattern = (
            _MBAT_START_PATTERN
            if self.method == "mbat_combo" else _FASTBAT_START_PATTERN
        )
        phases = [
            (
                "load_variants",
                None,
                "Load and filter LD-reference and GWAS variants",
            ),
            (
                "match_variants",
                _VARIANT_MATCH_PATTERN,
                "Match GWAS variants to the LD reference",
            ),
        ]
        if self.method in {"fastbat_gene", "mbat_combo"}:
            phases.append((
                "map_genes",
                _GENE_MAPPING_PATTERN,
                "Map gene positions to LD-reference variants",
            ))
        phases.append((
            "analyse_units",
            analysis_pattern,
            analysis_titles[self.method],
        ))
        self._phases = tuple(phases)
        self._phase_index = 0

    @property
    def enabled(self) -> bool:
        return self._progress.enabled

    @property
    def title(self) -> str:
        return self._phases[self._phase_index][2]

    def start(self) -> None:
        self._progress.start(self.title, total=self.total)
        self._record_phase()

    def _record_phase(self) -> None:
        phase, _, title = self._phases[self._phase_index]
        self.logger.record(
            "STATUS",
            "gcta_execution_phase",
            method=self.method,
            phase=phase,
            title=title,
            status="RUNNING",
        )

    def _advance_phases(self, text: str) -> None:
        for index in range(self._phase_index + 1, len(self._phases)):
            _, pattern, title = self._phases[index]
            if pattern is None or pattern.search(text) is None:
                continue
            self._phase_index = index
            self._progress.set_phase(title)
            self._record_phase()

    def _read_native_progress(self) -> None:
        text = self._native_log.read()
        self._advance_phases(text)
        if self._uses_result_rows:
            if self.total is None:
                for pattern in _MBAT_TOTAL_PATTERNS:
                    match = pattern.search(text)
                    if match is not None:
                        self.total = int(match.group(1))
                        break
            assert self._results is not None
            self.completed = self._results.count()
            return

        matches = _FASTBAT_PROGRESS_PATTERN.findall(text)
        if not matches:
            return
        native_completed, native_total = (int(value) for value in matches[-1])
        if self.total is None:
            self.total = native_total
        elif self.total != native_total and not self._native_total_mismatch_recorded:
            self.logger.record(
                "WARNING",
                "gcta_progress_native_total_mismatch",
                method=self.method,
                configured_work_units=self.total,
                native_work_units=native_total,
            )
            self._native_total_mismatch_recorded = True
        self.completed = max(
            self.completed,
            min(native_completed, self.total),
        )

    def refresh(self) -> None:
        self._read_native_progress()
        if self.total is not None and not self._total_recorded:
            self.logger.record(
                "OBSERVED",
                "gcta_progress_plan",
                method=self.method,
                units=self.unit_label,
                total=self.total,
                source=self._total_source,
            )
            self._total_recorded = True
        displayed = self.completed
        if self.total is not None:
            displayed = min(displayed, self.total)
        self._progress.update(displayed, total=self.total, title=self.title)
        if self.total is None:
            return
        safe_completed = min(displayed, max(0, self.total - 1))
        percentage = int(100 * safe_completed / self.total)
        if percentage == self._last_logged_percentage:
            return
        self._last_logged_percentage = percentage
        self.logger.record(
            "STATUS",
            "gcta_unit_progress",
            method=self.method,
            units=self.unit_label,
            completed=safe_completed,
            total=self.total,
            percentage=percentage,
            status="RUNNING",
        )

    def complete(self, validated_units: int) -> None:
        validated_units = int(validated_units)
        if validated_units < 1:
            raise ValueError("Validated GCTA progress requires at least one result unit")
        self.refresh()
        planned = self.total or validated_units
        if planned != validated_units:
            self.logger.record(
                "WARNING",
                "gcta_progress_result_count_difference",
                method=self.method,
                planned_work_units=planned,
                validated_result_units=validated_units,
            )
        self.completed = planned
        self.total = planned
        result_label = "result row" if validated_units == 1 else "result rows"
        self._progress.complete(
            planned,
            total=planned,
            title=(
                "Validate %s %s for %s"
                % (f"{validated_units:,}", result_label, self.unit_label)
            ),
        )
        self.logger.record(
            "STATUS",
            "gcta_unit_progress",
            method=self.method,
            units=self.unit_label,
            completed=planned,
            total=planned,
            validated_result_units=validated_units,
            percentage=100,
            status="VALIDATED",
        )

    def fail(self) -> None:
        self.refresh()
        self._progress.fail(title=self.title)
        completed = self.completed
        if self.total is not None:
            completed = min(completed, max(0, self.total - 1))
        self.logger.record(
            "STATUS",
            "gcta_unit_progress",
            method=self.method,
            units=self.unit_label,
            completed=completed,
            total=self.total,
            status="FAILED",
        )


def gcta_native_log_path(output_prefix: Path) -> Path:
    """Return GCTA's documented ``--out`` prefix plus ``.log`` path."""
    return Path(str(output_prefix) + ".log")


__all__ = ["GctaResultProgress", "gcta_native_log_path"]
