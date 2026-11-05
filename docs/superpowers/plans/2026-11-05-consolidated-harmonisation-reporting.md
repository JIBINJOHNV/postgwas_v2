# Consolidated Harmonisation Reporting Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [x]`) syntax for tracking.

**Goal:** Report input validation and completed study-wide decisions once, with one pooled genome-build decision and no routine chromosome build listing.

**Architecture:** Keep existing analysis/evidence producers and shared presentation. Consolidate input display, then move the detailed study report to the end of preparation. Existing failure reporting renders the available evidence once at a failing boundary; routine automatic inference notices move into that report only for orchestrated runs.

**Tech Stack:** Existing Python >=3.10, Polars, pytest, shared `postgwas.core.ui`, screen recorder and validation/file-reference presentation. Pinned local validation runtime: Python 3.10, Polars 0.20.31.

**Spec:** ../specs/2026-11-05-consolidated-harmonisation-reporting-design.md

## Global Constraints

- Read owning AGENTS.md and current cumulative diff before edits; preserve the completed AF implementation at `cc2bd07` and unrelated user files.
- No scientific calculations, thresholds, filtering, source-selection, AF approval contracts, stage order or checkpoint changes.
- Preserve internal chromosome build evidence and existing fallback-only veto semantics; remove routine per-chromosome build presentation only.
- Preserve input versus conversion/inference denominators, supplied external sources, warning/failure severity and known partial evidence.
- Keep mandatory progress and final output/completeness reporting; preparation success does not establish final analysis success.
- Reuse current renderers/logger/recorder/file references. No text-based/global deduplication or additional reporting framework. Keep presentation state outside scientific decisions and approval identities.
- No new scientific YAML keys/defaults, dependencies or edits to protected adapters.
- Implementation and scoped local verification are complete; exact server-scale validation remains unperformed.

## Review Focus

- Late external-AF or replacement failure after build/effect/P-value inference must yield one report with completed earlier evidence (Task 3).
- INFO inference, partition-write and final chromosome-set errors must also yield measured partial evidence once (Task 3).
- A normal build decision with recorded chromosome conflict evidence must retain its existing scientific outcome; hiding rows must not change fallback gates (Task 2).
- Hidden-screen recording and two datasets with identical decisions must retain one independent report per dataset, with no state leakage (Task 4).
- Different input and conversion denominators or duplicate filenames must remain distinguishable after consolidation (Tasks 1 and 4).

## File map

Paths below are relative to `/Users/JJOHN41/Downloads/gwas_sumstat/.implementation/postgwas_v2`.

| File | Responsibility |
|---|---|
| `src/postgwas/modules/harmonisation/service.py` | Input block composition and final/failure report ownership |
| `src/postgwas/modules/harmonisation/study_reporting.py` | Single build/frequency/INFO/inference presentation from saved evidence |
| `src/postgwas/modules/harmonisation/study_properties.py` | Orchestration-only deferral of routine auto-inference notices |
| `src/postgwas/modules/harmonisation/genome_build.py` | Stop routine readable chromosome observations; preserve full returned/extra evidence |
| Existing shared UI/logging/recorder | Reuse unchanged unless tests expose a required presentation-routing correction |
| Existing reporting, completeness, warning, build and transcript tests | Regression coverage described below |
| `docs/wiki/harmonisation/processing-order.md`, `outputs-and-qc.md` | Describe report timing and retained structured evidence |

### Task 1: One input-validation section

**Files:** `service.py` existing `_field_completeness_block`, `_input_validation_summary_block`, `_validated_column_mapping_block`; `tests/test_harmonisation_completeness.py`, `test_harmonisation_study_report_integration.py`.

**Interfaces:** Extend `_field_completeness_block(report, *, validation_summary=None, sample_columns=None, study_decisions=None)` to accept already-measured metadata/ledger values. Existing direct one-argument callers remain supported. Reuse the ledger formatter without printing a second top-level block. Restrict the later sample-sheet display to information not already shown, or retain it only as structured resolved configuration.

- [x] Write failing tests for one input section containing source/column mapping, field-level missing/conversion counts and one removal ledger. Use the pasted 11,749,765 input / 11,434,316 retained / 315,449 removed counts in pure reporting fixtures; assert they retain their labels/scopes and are not summed twice.
- [x] Add controls for user-excluded chromosomes, unequal physical/parsed row counts, configured-but-absent columns, external AF/INFO, unavailable Z and sample-count rounding. Two sources with the same basename must retain distinct established file references; unregistered templates retain their exact declaration once without new file reads.
- [x] Run the new tests RED; consolidate the existing renderers/call sites without rescanning input or recomputing counts. Use established file numbers and literal column mappings.
- [x] Verify counts/source identities and existing completeness artifacts remain unchanged; run focused tests GREEN and review the cumulative diff before committing this task.

### Task 2: One dataset build view and no repeated AF statistics

**Files:** `study_reporting.py`, `genome_build.py`; `tests/test_harmonisation_study_reporting.py`, `test_harmonisation_build_fallback.py`, `test_harmonisation_effective_frequency_reporting.py`.

**Interfaces:** Keep `infer_genome_build()` return and `step.extra['genome_build']` unchanged. `_build_fields(info, references)` displays pooled fields only. `_frequency_fields()` and the existing effective-screen formatter assign one owner to each saved numerator/denominator/cutoff; retain distinct missing/range and comparison evidence.

- [x] Add failing renderer tests with 22 saved chromosome records: one selected-build field, no `Chromosome N build evidence` rows, unchanged `info['per_chromosome']` and global metrics. Test ambiguous, declared, normal and fallback paths, including a normal-path fixture with recorded chromosome veto evidence.
- [x] Add a canonical-log check that routine `Genome-build chromosome N` observations are absent while full chromosome evidence remains in stage extras and successful/failed manifests.
- [x] Add a native-AF test proving each usable/unusable count, low-frequency fraction and cutoff is shown once, with A and B independently labelled when both exist.
- [x] Run RED, remove the chromosome display/observation loops and duplicate AF display fields, then run build/science/reporting controls GREEN. Do not remove chromosome evidence computation or change gate outcomes. Review cumulative diff and commit.

### Task 3: One final preparation report, including failures and warnings

**Files:** `service.py`, `study_properties.py`, `study_reporting.py`; `tests/test_harmonisation_study_report_integration.py`, `test_harmonisation_replacement_af_reporting.py`, `test_harmonisation_pvalue_inference_warning.py`, `test_harmonisation_effect_detection.py`.

**Interfaces:** Preserve `render_study_decisions(...) -> str`, `_study_decisions_block(...)` and `_report_study_failure(...)` as shared evidence-based presentation. Add keyword-only `report_inference_warnings=True` to `resolve_study_properties()`; orchestration passes `False`. Thread this only to routine automatic effect/P-value notices; keep saved warning strings and direct defaults. Fatal/declaration-mismatch warnings remain immediate. Extend the consolidated report with saved INFO type/source/evidence, using existing field helpers.

- [x] Add failing real-preparation tests for success, external-AF failure after partitioning, replacement rejection, staging/source failure, INFO inference failure, partition-write failure and final chromosome-set mismatch. Assert exactly one expanded study report and one known build decision in each run; no pending report followed by a second full report.
- [x] Add tests that routine effect/P-value warning evidence, reason and recommendation appear once in the orchestrated report, direct callers still receive notices, and declaration mismatch/fatal diagnostics remain immediate and visible on failure.
- [x] Run RED. Remove the pre-step-8 full render and standalone AF/replacement cards from orchestration. Render the consolidated report once after AF/INFO preparation and chromosome reconciliation, immediately before returning prepared worker inputs. Add missing failure-boundary coverage using the existing failure helper, preserving original exception types/partial evidence.
- [x] For a failure after this report has already been emitted, retain the failing-stage/error summary without replaying build/effect/P-value sections. Avoid text matching and scientific-state flags; ownership belongs to the preparation lifecycle.
- [x] Update the existing provisional-output ordering assertion intentionally: provisional A/B evidence remains recorded internally, but only the terminal preparation report appears in readable output. Keep all numeric gate/provenance assertions.
- [x] Run integration/warning/report suites GREEN; verify frames, AF contracts, selected sources and genome-build evidence are unchanged. Review cumulative diff and commit.

### Task 4: Transcript regression, documentation and handoff

**Files:** `tests/test_harmonisation_study_reporting_transcripts.py`, `test_harmonisation_logging.py`, existing AF/build/completeness regressions; the two wiki pages in the file map.

**Interfaces:** Exercise real orchestration through existing `record_screen`, `_DatasetScreenRouter` and `PipelineLogger`, rather than only manually calling a renderer. Keep native diagnostics and shared progress behavior intact.

- [x] Add visible/hidden and redirected transcript tests for one and two datasets, success and late failure, narrow/wide wrapping and active progress. Assert one input section, one expanded study section, one selected build, zero routine chromosome build rows and no standalone repeated inference cards in each relevant log surface.
- [x] Assert early failures have known partial evidence only; failures never show unperformed AF/INFO approval. Preserve warning severity and literal error/path text. Verify renderer calls never open inputs or invoke build/AF inference/joins.
- [x] Run focused reporting/build/completeness/warning tests with `POLARS_MAX_THREADS=4 PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src /Users/JJOHN41/miniconda3/envs/postgwas/bin/python -B -m pytest -q -p no:cacheprovider` plus the file-map test paths. Verify numerical/source contracts against current AF regression suites, including `test_harmonisation_unified_frequency_policy.py`.
- [x] Update the processing/output docs; run policy/wiki checks and the broader selected regression suite. Report existing unrelated collection/historical-baseline/CLI-documentation issues separately; do not hide or weaken them.
- [x] Save representative before/after transcripts and exact command/results. Review the complete cumulative diff, apply AGENTS.md compliance checks, preserve the original source data and remove task-created temporary fixtures before handoff.

## Planned output order

`Dataset identity / consolidated input validation` → stage progress → `Study-wide decisions and evidence` (build, effect, SE, P-value, BETA/Z, strand, AF A/B, INFO) → chromosome progress and actual row QC → final output/completeness results.

If preparation fails, its one study report ends with the failing check and dataset stop reason. A stage-5 start/completion line remains progress, not a second build decision.

## Plan verification

Execution completed with 1,394 selected tests passing across 49 files in three bounded batches (546 + 588 + 260). Independent whole-patch review found no blocking issue; 61 focused checks verified preserved canonical warning counts and failure/transcript handling. Policy reference and 50-page wiki checks passed. Existing unrelated collection/historical-summary/CLI-extraction issues and exact server-scale validation are documented in `/Users/JJOHN41/Downloads/gwas_sumstat/harmonisation_reporting_implementation_20261105/validation.md`.

Publication was explicitly authorized by the user. Tasks share presentation interfaces and are saved as one atomic reporting commit after verification; the separately verified cc2bd07 AF commit is its prerequisite. Scientific configuration, source selection, thresholds and protected adapters remain unchanged.
