# Dual whole-study AF rules implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [x]`) syntax for tracking.

**Goal:** Implement and commit the approved two-rule whole-study AF decision.

**Architecture:** Extend the existing pooled comparison and owning study preflight;
keep prepared-match reuse and all ordinary external-source behavior. Isolate
whole-study thresholds in canonical YAML and expose both rules through existing
evidence/reporting interfaces. One integrated task keeps source, config, tests
and documentation consistent; independent task and final reviews follow.

**Tech Stack:** Python, Polars >=0.20, canonical YAML policy registry, Pydantic,
pytest. Use the installed postgwas conda Python; do not install dependencies.

**Spec:** [Approved design](../specs/2026-11-04-dual-whole-study-af-rules-design.md)

## Global constraints

- Follow AGENTS.md before and after every change cycle; preserve unrelated files.
- Preserve existing MAF-like entry screen and source rows/prepared matches.
- Require >=1,000 usable pairs and >=0.50 coverage, with a fixed pre-exclusion denominator.
- Both full rules run automatically after coverage passes on the same pairs.
- Rule 1: corr(x,r)>=0.95; mean(abs(x-r))<=0.10; P(r<=0.5)>0.90.
- Rule 2: corr(x,min(r,1-r))>=0.95; mean(abs(x-min(r,1-r)))<=0.10; P(r>0.5)>0.20.
- Full-rule truth table: TF retain; FT replacement; TT and FF inconclusive/stop.
- No whole-study 0.02 error gap, alternative ambiguity exception, retry or informative-subset classifier.
- Replacement requires both resolved public settings; use shared required-argument validation.
- Canonical YAML is the single default source. Keep ordinary external thresholds/behavior unchanged.
- No protected adapters, new dependencies, source-data mutation or extra reference reads.
- Document heuristics and uncertainty; do not claim production-scale calibration.

## Review focus

- Raw supplied AF versus folded reference MAF must not accidentally fold source AF.
- Coverage failure must prevent either interpretation/replacement from authorizing continuation.
- Both-pass guard needs a real configurable overlapping-profile control, not default unreachable assumptions.
- Reports, worker decisions, provenance and checkpoint identities must agree on both-rule outcomes.
- Ordinary external classification and non-MAF-like study inputs must retain their previous behavior.

## Task 1: Implement both complete rules and their integrated interfaces

**Files:**
- Modify: src/postgwas/modules/harmonisation/allele_frequency.py
- Modify: src/postgwas/modules/harmonisation/study_reporting.py
- Modify: src/postgwas/config/defaults/modules/harmonisation.yaml
- Modify if necessary: src/postgwas/config/models/modules/harmonisation.py
- Modify: tests/test_harmonisation_af_acceptance.py and affected tracked config/report tests
- Modify: docs/modules/harmonisation/configuration.md
- Modify: docs/modules/harmonisation/eaf_replacement.md
- Modify: docs/wiki/harmonisation/policies.md

**Interfaces:** Extend _compare_frequency_pairs and its existing whole-study caller
without changing ordinary callers. Whole-study evidence must expose both rule
booleans, raw supplied-to-reference-MAF correlation, reference-major fraction,
both errors, failed requirements and a final combined reason. Keep legacy evidence
meaning clear; the existing maf_correlation denotes folded compatibility.
Use canonical new eaf.study_reference_min_correlation=0.95 and
eaf.study_reference_min_major_fraction=0.20; register in module/reporting and
preparation policy identities. Keep eaf.maf_reference_min_correlation=0.90 and
eaf.maf_reference_error_margin=0.02 for ordinary external callers.

- [x] Read every affected implementation, caller, config, test and consumer before editing.
- [x] Write meaningful failing tests first; record RED commands/results in the task report.
  Test actual whole-study retain, MAF/replacement, both-fail and both-pass cases;
  use configurable overlapping profile bounds for the both-pass case. Verify
  rule 2 is measured even on rule 1 pass. Test source [0.2,0.8]-like controls
  to distinguish raw-source correlation from folded-source compatibility.
- [x] Test exact count/coverage/correlation/error/profile boundaries, empty/missing/
  invalid/excluded pairs, constant direct and MAF vectors, and the 10–20% gap.
  Cover complete-failure reasons and no retry/replacement for coverage failure.
- [x] Test both missing replacement settings, each missing partner, invalid source
  validation, original AF/provenance preservation and chromosome reuse/orientation.
- [x] Test canonical defaults/export/reload/invalid settings and ordinary external
  classifier unchanged. Valid numeric-string policy overrides remain supported.
- [x] Implement minimal extensions through existing helpers and single aggregation.
  No new .get fallback thresholds. Keep the final combined decision as the sole
  authority for accept/replacement. Preserve the agreed branch on partial EAF
  agreement plus full MAF pass; no added exception.
- [x] Update reporting and all directly affected consumers, evidence and docs.
  Describe global QC support, not semantic proof; keep the documented known
  near-half/population-difference counterexample as a limitation, not a new guard.
- [x] Run focused tests GREEN and inspect complete cumulative diff. Run shared
  config/export and affected reporting/checkpoint/interface regressions. Update
  superseded expected behavior only where the approved spec intentionally changes it.
- [x] Controller runs every Git-tracked test and bare pytest; report existing
  unrelated collection errors explicitly rather than widening scope to fix them.
- [x] Review hardcoding/duplication, source/protected-file safety and generated
  artifacts; report AGENTS compliance fields. Do not claim empirical calibration.
- [x] Controller obtains independent task and whole-change reviews and resolves findings.

Delivery: stage only authorized files and create the user-requested local commit.
No merge or push is requested.

## Validation record

Implementation follows the approved user-selected QC rules. Their exact cutoffs
have not been empirically calibrated; these tests validate implementation and
numerical behavior, not universal identification of frequency provenance.

### Environment and commands

Python: 3.10.21
Polars: 0.20.31
Pydantic: 2.13.4
PyYAML: 6.0.3
pytest: 8.4.2

Every pytest command used the installed postgwas conda Python with
`PYTHONDONTWRITEBYTECODE=1`, `PYTHONPATH=src`, the conda bin directory first in
`PATH`, `python -B -m pytest -p no:cacheprovider`, and an owned temporary directory.
No new dependencies were installed.

| Verification | Expected and observed result |
|---|---|
| New regression tests before implementation | RED: 34 failed / 89 passed; separate exact-boundary RED: 7 failed. Missing new policy/evidence and absent branches caused the expected failures. |
| Focused release: AF acceptance, module config, config export, reference reuse | GREEN: 240 passed in 115.00s after final logging/help fixes. |
| Extended AF/reporting/replacement/config checks | 274 passed / 1 old Spearman expectation failed in 150.63s. The fixture's 0.928609 correctly fails the new 0.95 default. Its revised test checks default rejection and explicit study-only 0.90 acceptance. Full real-record module then passed: 32 tests in 9.28s. |
| Earlier all 14 Git-tracked test files | 512 passed, 42 existing Polars join deprecation warnings, 216.71s. |
| Final scientific source: all 14 Git-tracked test files | **513 passed**, 42 existing Polars join deprecation warnings, **225.41s**. |
| Final two-row HTML reporting follow-up | Expected RED: 3 failures. GREEN: 4 HTML regressions in 0.79s and 36 affected reporting/integration/transcript tests in 14.07s. |
| Public config export and reload | New whole-study defaults 0.95 / 0.20 exported and validated. Ordinary external defaults 0.90 / 0.02 preserved. |
| Public config invalid values | Both invalid new policy keys reported together; exit 2; no traceback. |
| Policy documentation generator | `python -B tools/docs/update_harmonisation_policies.py --check`: synchronized. |
| Source/configuration hygiene | Changed Python ASTs parse; YAML loads; `git diff --check` passes; no protected adapter diff or new dependencies. No formatter/linter/type-check command is configured in pyproject.toml. |
| Bare project pytest | **Blocked by 37 pre-existing collection errors**, 1 existing marshmallow/distutils warning, 9.66s; details below. |

The final tracked-suite command selected the exact paths returned by
`git ls-files tests` and invoked pytest with those 14 paths. The unrestricted
command was the same pytest prefix without selected test paths. Focused release
paths were `tests/test_harmonisation_af_acceptance.py`,
`tests/test_harmonisation_config.py`, `tests/test_harmonisation_config_export.py`
and `tests/test_harmonisation_frequency_reference_reuse.py`.

The tests cover all four combined complete-rule outcomes; actual replacement
selection and both-pass rejection; raw supplied versus folded reference
correlation; strict 10–20% profile gap; exact numerical boundaries; insufficient,
empty, excluded and constant evidence; required/invalid replacement settings;
original-value preservation; worker reuse and single swaps; canonical failures,
terminal/HTML reporting; config exports; and policy-dependent checkpoint and
preparation identities. Targeted RED/GREEN also demonstrated and corrected
supported MAF being labelled EAF on a missing-setting failure, duplicate fatal
messages, duplicate policy logging, missing HTML evidence and stale preparation
under whole-study policy changes.

One initial baseline test under the Documents checkout's temporary directory
reported a changed reference identity on chromosome 2 (143 other tests passed).
It did not recur in an independent diagnostic run of all 84 original AF tests
under private temporary storage, or either complete tracked regression run.
The diagnostic compared every contract section and observed no drift. Strict
file fingerprints were preserved; no reproducible runtime defect was found.

### Repository scope and implementation decisions

Changed runtime/configuration files: `allele_frequency.py`, `strand.py`,
`study_reporting.py`, `html_report.py` and canonical harmonisation YAML. Updated
the tracked AF acceptance regression file, module configuration/replacement
reference documentation, generated policy reference, this plan and its spec.
The new canonical controls are `eaf.study_reference_min_correlation=0.95` and
`eaf.study_reference_min_major_fraction=0.20`. Legacy external controls remain
available with their unchanged defaults.

Two already ignored local historical test files have narrow expectation updates
for the approved semantics: `tests/test_harmonisation_af_failure_reporting.py`
and `tests/test_harmonisation_study_eaf.py`. They remain local under the existing
repository ignore policy; corresponding self-contained regression coverage is
in the tracked AF acceptance test. Existing `_to_delete/` and `tmp/` directories
are unrelated and preserved.

Rulings made during implementation:

1. Use a distinct whole-study correlation key to preserve ordinary external
   behavior. If incorrectly scoped, an existing external classifier could
   otherwise receive an unintended threshold change.
2. Reuse the clean existing feature checkout. No worktree, branch switch or
   unrelated-file cleanup is needed; the risk of an incorrect choice would be
   overlapping user edits, which the cumulative diff audits checked.
3. Proceed under the user's explicit approval of the chat plan and local commit;
   do not request the same approval again. Routine implementation remains
   reviewable in the local diff and commit.
4. Bind every whole-study AF decision policy into the prepared-match identity,
   not just the two new keys. This prevents stale acceptance under changed
   library policies; conservative invalidation can require a safe restart.

### Full-project collection blockers

These paths were not changed by this task. Thirty-six errors are duplicate
module names in copied output snapshots and the original tools tree. The
coordinate test imports `_stage_indexed_csv_batches`, which does not exist in
its unchanged concordance-service dependency. No collection rule was weakened
and no unrelated snapshot or test was deleted.

- `output/installation_validation_20261030/source_candidate_02/tools/resource_preparation/dbsnp_alfa_latest/test_prepare_population_frequencies.py`
- `output/installation_validation_20261030/source_candidate_02/tools/resource_preparation/dbsnp_build157_grch37/test_prepare.py`
- `output/installation_validation_20261030/source_candidate_02/tools/resource_preparation/finngen_r13_annotations/test_prepare_info_scores.py`
- `output/installation_validation_20261030/source_candidate_02/tools/resource_preparation/thousand_genomes_phase3/test_prepare.py`
- `output/installation_validation_20261030/source_candidate_03/tools/resource_preparation/dbsnp_alfa_latest/test_prepare_population_frequencies.py`
- `output/installation_validation_20261030/source_candidate_03/tools/resource_preparation/dbsnp_build157_grch37/test_prepare.py`
- `output/installation_validation_20261030/source_candidate_03/tools/resource_preparation/finngen_r13_annotations/test_prepare_info_scores.py`
- `output/installation_validation_20261030/source_candidate_03/tools/resource_preparation/thousand_genomes_phase3/test_prepare.py`
- `output/installation_validation_20261030/source_candidate_04/tools/resource_preparation/dbsnp_alfa_latest/test_prepare_population_frequencies.py`
- `output/installation_validation_20261030/source_candidate_04/tools/resource_preparation/dbsnp_build157_grch37/test_prepare.py`
- `output/installation_validation_20261030/source_candidate_04/tools/resource_preparation/finngen_r13_annotations/test_prepare_info_scores.py`
- `output/installation_validation_20261030/source_candidate_04/tools/resource_preparation/thousand_genomes_phase3/test_prepare.py`
- `output/installation_validation_20261030/source_candidate_05/tools/resource_preparation/dbsnp_alfa_latest/test_prepare_population_frequencies.py`
- `output/installation_validation_20261030/source_candidate_05/tools/resource_preparation/dbsnp_build157_grch37/test_prepare.py`
- `output/installation_validation_20261030/source_candidate_05/tools/resource_preparation/finngen_r13_annotations/test_prepare_info_scores.py`
- `output/installation_validation_20261030/source_candidate_05/tools/resource_preparation/thousand_genomes_phase3/test_prepare.py`
- `output/installation_validation_20261030/source_candidate_06/tools/resource_preparation/dbsnp_alfa_latest/test_prepare_population_frequencies.py`
- `output/installation_validation_20261030/source_candidate_06/tools/resource_preparation/dbsnp_build157_grch37/test_prepare.py`
- `output/installation_validation_20261030/source_candidate_06/tools/resource_preparation/finngen_r13_annotations/test_prepare_info_scores.py`
- `output/installation_validation_20261030/source_candidate_06/tools/resource_preparation/thousand_genomes_phase3/test_prepare.py`
- `output/installation_validation_20261030/source_candidate_07/tools/resource_preparation/dbsnp_alfa_latest/test_prepare_population_frequencies.py`
- `output/installation_validation_20261030/source_candidate_07/tools/resource_preparation/dbsnp_build157_grch37/test_prepare.py`
- `output/installation_validation_20261030/source_candidate_07/tools/resource_preparation/finngen_r13_annotations/test_prepare_info_scores.py`
- `output/installation_validation_20261030/source_candidate_07/tools/resource_preparation/thousand_genomes_phase3/test_prepare.py`
- `output/installation_validation_20261030/source_initial/tools/resource_preparation/dbsnp_alfa_latest/test_prepare_population_frequencies.py`
- `output/installation_validation_20261030/source_initial/tools/resource_preparation/dbsnp_build157_grch37/test_prepare.py`
- `output/installation_validation_20261030/source_initial/tools/resource_preparation/finngen_r13_annotations/test_prepare_info_scores.py`
- `output/installation_validation_20261030/source_initial/tools/resource_preparation/thousand_genomes_phase3/test_prepare.py`
- `output/sample_sheet_metadata_validation_20261030/source/tools/resource_preparation/dbsnp_alfa_latest/test_prepare_population_frequencies.py`
- `output/sample_sheet_metadata_validation_20261030/source/tools/resource_preparation/dbsnp_build157_grch37/test_prepare.py`
- `output/sample_sheet_metadata_validation_20261030/source/tools/resource_preparation/finngen_r13_annotations/test_prepare_info_scores.py`
- `output/sample_sheet_metadata_validation_20261030/source/tools/resource_preparation/thousand_genomes_phase3/test_prepare.py`
- `tests/test_harmonisation_coordinates.py`
- `tools/resource_preparation/dbsnp_alfa_latest/test_prepare_population_frequencies.py`
- `tools/resource_preparation/dbsnp_build157_grch37/test_prepare.py`
- `tools/resource_preparation/finngen_r13_annotations/test_prepare_info_scores.py`
- `tools/resource_preparation/thousand_genomes_phase3/test_prepare.py`

### Limits and final review

No full production GWAS rerun, external-tool/server integration run, threshold
calibration or 100–200-million-row benchmark was performed. Existing Polars
coalesce warnings come from unchanged joins; full project collection remains
blocked as above. Focused and tracked-suite results do not establish empirical
classification accuracy or production-scale memory/runtime.

Independent task review approved spec compliance and code quality with no
Critical or Important findings. Final whole-change review identified one minor
HTML omission: failed shared coverage gates hid observed pair counts and the
fixed denominator because profile metrics were unavailable. A two-row
presentation-only correction now renders the recorded counts independently of
those metrics. Three real failed-comparator/full-renderer regressions first
failed as expected, then all four HTML cases and 36 reporting checks passed.
The scientific algorithm, parameters and source evidence were unchanged by this
follow-up, so the 513-test broad suite was not repeated for those two HTML rows.
Scoped final re-review marked the finding addressed and found no new breakage.

Checklist: passed for the scoped implementation and reporting corrections.
Cumulative diff reviewed: yes. Hardcoding and duplication audit: passed.
Scientific validation: definitions, deterministic numerical controls and
representative fixture behavior checked; empirical calibration and production
scale remain unverified. Tests and exclusions are documented above. Outstanding
repository-wide issues are the 37 pre-existing collection errors and existing
warnings; no unresolved finding remains in the authorized patch. This validated
record supports the user-authorized local commit. No merge or push is requested.
