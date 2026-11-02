# Study-wide report implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox syntax for tracking.

**Goal:** Display the evidence behind study-wide decisions together without repeating scientific analysis.

**Architecture:** A pure harmonisation renderer consumes existing decision dictionaries and resource identities. Service orchestration renders resolved decisions after AF validation and renders partial evidence on failure; existing calculations and failure policy stay authoritative.

**Tech Stack:** Python, existing Polars calculations, shared PostGWAS terminal/logging helpers, pytest.

**Spec:** `docs/superpowers/specs/2026-11-02-study-wide-report-design.md`.

## Global constraints

- Work on main; preserve pending compressed-input changes and unrelated files.
- Protected adapters, scientific algorithms, thresholds and configuration defaults remain unchanged.
- Reuse existing scalar evidence; no extra input/reference reads for reporting.
- Tests remain local under the repository's existing ignored-tests policy.
- Test with available original real data before pushing; do not claim a full native run from a bounded validation.

## Review focus

- No-comparison frequency acceptance must not claim reference confirmation.
- MAF-supported replacement selection must not claim completed EAF assignment.
- Missing inputs and empty eligible subsets must not produce fabricated passes/counts.
- Explicit strand declarations must not appear as inferred consensus.
- Early and AF failures must retain their actual evidence and batch behavior.

## Task 1: Pure evidence renderer and failure handoff

**Files:** Create `src/postgwas/modules/harmonisation/study_reporting.py`; modify the existing AF exception boundary in `allele_frequency.py`; add focused local renderer/failure tests.

**Interfaces:** `render_study_decisions(study_decisions, *, genome_build_info=None, build_reference_files=None, policies=None, resource_maps=None) -> str`; optional `StudyAlleleFrequencyError.study_decisions` contains copied existing evidence.

- [ ] Write and observe failing tests for conditional/missing/failed evidence and width alignment.
- [ ] Render using shared UI and existing scalar formatters; preserve scientific decision dictionaries.
- [ ] Attach already-computed AF failure evidence without altering classification or error messages.
- [ ] Run the focused tests and inspect the complete cumulative diff.

## Task 2: Orchestration, reproducible evidence and delivery

**Files:** Modify `service.py`, `study_properties.py`, focused local integration tests and the existing harmonisation documentation.

- [ ] Write and observe failing integration tests for resolved-report ordering, manifest evidence and early failures.
- [ ] Replace the short report boundary with the pure renderer, after conditional AF validation; retain partial study-property evidence on errors.
- [ ] Persist existing scalar evidence and keep progress, warnings and scientific outputs unchanged.
- [ ] Run affected tests, required documentation checks and the bare pytest collection attempt; report independent baseline failures explicitly.
- [ ] Run original real data through the affected scientific stages against a frozen pre-change package; compare decisions, row accounting and available data outputs.
- [ ] Review the complete diff, remove owned temporary artifacts, commit completed changes on main and push origin after the required validation.
