# Whole-study Replacement AF Validation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Validate a provisionally selected replacement AF column once for the
whole study, using the existing frequency attachment and the five agreed gates.

**Architecture:** Preserve the original dual-rule decision. Reuse the shared
external merge, staged reference projections and default orientation evidence;
add one pooled replacement gate before chromosome workers. Generalize existing
metric aggregation instead of invoking either existing classifier as a new rule.

**Tech Stack:** Existing Python >=3.10, Polars, PyArrow, canonical YAML/configuration,
shared logging/checkpoints and pytest; no new dependency is proposed.

**Spec:** `docs/superpowers/specs/2026-11-04-replacement-af-validation-design.md`

**Status (2026-11-04):** implemented and independently reviewed; final scoped
checks passed (782 workspace cases and 317 clean-candidate cases). The local-only
commit record and explicit baseline/scale limitations are below.

## Global Constraints

- Implement the approved design, preserving unrelated working-tree changes.
- Original A is unchanged; B always refers to the original listed study EA.
- Gates: n >=1,000; n/N >=0.50; corr(A,min(B,1-B)) >=0.95;
  MAE <=0.10; fraction(B>0.5) >0.20, using the spec's common paired set.
- Reuse existing canonical policies/method; no .02/.05 error-margin gate,
  additional reference, duplicated merge, new default or silent fallback.
- Retain existing multiallelic attachment and biallelic diagnostic eligibility.
- Preserve source row identity/order, original MAF, BETA/Z, default-reference
  orientation and original AF-decision evidence. Do not edit protected adapters.
- Preserve unrelated working-tree changes and apply the repository change-cycle
  gate before and after each implementation task.

## Review Focus

1. Raw reverse-complement or swapped study keys must not reverse the major profile.
2. Same-source reuse must check the entire contract and never claim independence.
3. Diagnostic exclusions must not shrink N or silently remove output variants.
4. Source selection must not allow workers before replacement validation passes.
5. Chunking, chromosome sizes and rank ties must not change global statistics.

---

### Task 1: Pooled replacement evidence and canonical policies

**Files:**
- Modify: `src/postgwas/modules/harmonisation/allele_frequency.py`
- Modify: `src/postgwas/config/defaults/modules/harmonisation.yaml`
- Modify: `tests/test_harmonisation_af_acceptance.py`
- Modify: `tests/test_harmonisation_eaf_replacement_science.py`

**Interfaces:**
- Consume raw A and original-EA-aligned B, fixed study row count N, resolved
  policies and the existing scientific step context.
- Extend the existing frequency-pair measurement implementation for reuse by
  original and replacement decisions; do not duplicate its correlation/MAE code.
- Produce a replacement evidence dictionary with `status`, method, source
  identity/mapping, n/N, excluded counts, thresholds, measured statistics and
  ordered failed requirements. `status` is `passed` or `inconclusive`; source
  and configuration errors retain their existing exceptions.

- [x] Read each owning implementation, configuration, caller and test completely;
  verify current defaults and cumulative diff before editing.
- [x] Add failing tests for pass/fail/equality of each gate, zero usable pairs,
  missing/non-finite/out-of-range values, constant vectors and all exclusions.
- [x] Assert raw A is not folded or mutated; A=.20 with B=.79 compares to .21,
  not .79. Test no legacy margin or original dual-rule truth table is invoked.
- [x] Run the new focused tests and verify their expected failures.
- [x] Generalize shared measurement/count handling and implement the new all-five
  decision. Reuse existing five policy keys and correlation-method key exactly
  as listed in the spec; extend canonical help, not argparse defaults.
- [x] Implement stable pooled Pearson statistics and test against the unchunked
  result. Preserve exact global Spearman ranks/ties under existing memory
  controls; never substitute chunk-local ranks or Pearson.
- [x] From the repository root in the analysis/test environment, run
  `PYTHONPATH=src python -m pytest -p no:cacheprovider tests/test_harmonisation_af_acceptance.py tests/test_harmonisation_eaf_replacement_science.py`;
  expected: all pass, including original dual-rule and ordinary external-classifier regressions.
- [x] Check exported/reloaded configuration, boundary overrides, invalid values
  and policy logging. Complete diff/compliance review before the next task.

### Task 2: Existing merge, staging, adoption and checkpoint integration

**Files:**
- Modify: `src/postgwas/modules/harmonisation/service.py`
- Modify: `src/postgwas/modules/harmonisation/allele_frequency.py`
- Extend only if needed: `src/postgwas/modules/harmonisation/strand.py`,
  `src/postgwas/modules/harmonisation/external_reference_staging.py`
- Test: `tests/test_harmonisation_eaf_replacement_resources.py`,
  `tests/test_harmonisation_eaf_replacement_downstream.py`,
  `tests/test_harmonisation_frequency_reference_reuse.py`,
  `tests/test_harmonisation_multiallelic_reference_frequency.py`

**Interfaces:**
- Add `validate_replacement_allele_frequency(chromosome_partitions,
  sample_columns, study_decisions, resource_maps, external_eaf_colmap, *,
  policies=None, logger=None) -> dict` in `allele_frequency.py`; return updated
  study decisions or raise `StudyAlleleFrequencyError` carrying the evidence.
- Consume typed partition paths and staged source maps after dataset step 08
  staging. Reuse validated original-reference orientation and the same shared
  merge. Store evidence at `eaf_replacement_validation`; leave
  `eaf_reference_validation` unchanged.
- Replacement adoption requires `eaf_replacement_validation.status == 'passed'`;
  the earlier `replacement_selected` state denotes provisional source selection.

- [x] Add failing integration tests: R1-only never opens candidate contents;
  R2-only selects provisionally, pass permits workers, gate failure permits no
  scientific worker/tool invocation, both/neither original rules stop as before.
- [x] Add tests for missing both/one settings, unreadable/malformed sources,
  wrong column/build/mapping and complete actionable canonical errors.
- [x] Verify direct, swapped, reverse-complement and combined orientations use
  original EA for B's major fraction and leave A/BETA/Z/source rows unchanged.
- [x] Test direct and swapped multiallelic attachment through the existing merge,
  complete-site sums, invalid/conflicting/incomplete components and diagnostics
  excluding multiallelic rows while keeping them in fixed N and source outputs.
- [x] Add equivalent whole-genome/chromosome-template tests; enforce one global
  outcome and existing chromosome attachment coverage (default >=80%) separately.
- [x] Add same-source alias/mapping/column reuse tests and changed source/policy
  checkpoint invalidation. Fully identical reuse must not reread/stage the source;
  partially reused shared sources retain the existing single staging pass needed
  by their other chromosomes.
- [x] Run new tests to verify expected failures, then integrate the gate inside
  dataset preparation after staging and before step completion/worker launch.
  Use narrow partition reads and retained orientation evidence; no full effect
  harmonisation rerun just to prepare diagnostic keys.
- [x] Retain or reuse aligned attachment results under exact row/source/policy
  contracts where possible; avoid duplicate source parsing and new whole-study
  wide frames. Never overwrite default-reference prepared candidates.
- [x] Route all failures through canonical study failure/progress finalization;
  partition/staging intermediates cannot be mistaken for completed results.
- [x] Run the four listed integration files plus replacement science/acceptance;
  expected: all pass. Review the full cumulative diff and complete compliance.

### Task 3: Reporting, documentation, release verification and scale evidence

**Files:**
- Modify: `src/postgwas/modules/harmonisation/study_reporting.py`,
  `src/postgwas/modules/harmonisation/html_report.py`
- Modify: `docs/modules/harmonisation/eaf_replacement.md`,
  `docs/modules/harmonisation/configuration.md`
- Regenerate only through its owner if affected:
  `docs/wiki/harmonisation/policies.md` using
  `tools/docs/update_harmonisation_policies.py`
- Test: `tests/test_harmonisation_af_acceptance.py` and the existing reporting,
  progress, configuration, resume and public CLI test owners found in pre-edit
  inspection; extend those owners, do not create parallel infrastructure.

**Interfaces:**
- Terminal/canonical/HTML consume the same replacement evidence dictionary.
- Report original AF decision and provisional/accepted replacement distinctly;
  source provenance remains population-reference-imputed EAF.

- [x] Add failing tests for all numeric gates, multiple failures, missing metrics,
  fixed denominator, source identities, hidden-screen logging and failure progress.
- [x] Implement report fields and document both decisions, exact denominators,
  config keys, multiallelic attachment versus diagnostic eligibility, and limits.
- [x] Record same-source evidence as reused; report that folding/profile agreement
  supports consistency under the declared allele mapping, not semantic proof.
- [x] Add a counterexample characterization: A=[.1,.2,.3,.4] repeated; B=[.1,.8,.3,.6]
  and 1-B both give corr=1, MAE=0 and major fraction=.5. This documents the
  mapping limitation without introducing an unapproved extra decision gate.
- [x] Verify canonical config export/reload and exact threshold messages; no new
  public flags are required. Regenerate affected owned documentation and inspect.
- [x] Run scoped tracked-plus-feature regression and clean-candidate tests,
  applicable syntax/configuration/documentation checks, and independent review.
  Record pre-existing full-suite and documentation-check blockers explicitly.
- [x] Run deterministic chunk/partition-size invariance tests for Pearson and
  global Spearman with ties, highly unequal chromosome sizes, changed row order
  and missing data. Expected global statistics match the unchunked oracle within
  a documented numerical tolerance and give the same boundary decisions.
- [x] Record numeric-only 100M/200M and smaller shared-merge benchmarks, wall
  time and sampled RSS, with explicit omissions: OS peak memory, production
  I/O/spill and full-workflow capacity were not measured. No production run
  or 100–200M full-workflow capacity claim is authorized by these measurements.
- [x] Review complete cumulative diff, scientific assumptions, all call sites,
  hardcoding/duplication, checkpoint dependencies and protected paths. Clean only
  task-created temporary artifacts and report every unresolved test/scale limit.
- [ ] Commit implementation only under the user's implementation authorization;
  report the commit, files/config changes, observed versus expected validation,
  tests not run and scientific limits. No push is implied.

## Final implementation and validation record (2026-11-04)

Implementation and independent review are complete; final scoped verification
passed. P1 (alias retargeting) and P2 (private fixture dependency), plus the
actively imported fixture-path export correction, are addressed and approved.
The authorized local commit is the remaining repository operation: stage only
the 21 paths below after inspecting the exact staged diff; no push is authorized.
No commit hash is invented in this pre-commit record. Git records the resulting
local commit. Unrelated exponentiated-SE work in `e943add` is preserved and is
not attributed to this feature.

### Final scope: 21 owned files

```text
docs/modules/harmonisation/configuration.md
docs/modules/harmonisation/eaf_replacement.md
docs/wiki/harmonisation/policies.md
src/postgwas/config/defaults/modules/harmonisation.yaml
src/postgwas/modules/harmonisation/allele_frequency.py
src/postgwas/modules/harmonisation/html_report.py
src/postgwas/modules/harmonisation/service.py
src/postgwas/modules/harmonisation/shared/runtime.py
src/postgwas/modules/harmonisation/strand.py
src/postgwas/modules/harmonisation/study_reporting.py
tests/test_harmonisation_af_acceptance.py
tests/test_harmonisation_frequency_reference_reuse.py
docs/superpowers/plans/2026-11-04-replacement-af-validation.md
docs/superpowers/specs/2026-11-04-replacement-af-validation-design.md
tests/test_harmonisation_replacement_af_validation.py
tests/test_harmonisation_replacement_af_integration.py
tests/test_harmonisation_replacement_af_reporting.py
tests/test_harmonisation_eaf_replacement_resources.py
tests/test_harmonisation_eaf_replacement_science.py
tests/test_harmonisation_batch_preflight.py
tests/test_harmonisation_resource_preflight.py
```

Canonical YAML changes are help/origin text only. No configuration key, value,
default, public flag or scientific cutoff was added or changed. The six reused
AF defaults were checked unchanged. No protected adapter was edited.

### Final environment and observed tests

Verified analysis environment: Python 3.10.21, Polars 0.20.31, NumPy 1.23.5,
bcftools/htslib 1.23.1; Darwin arm64, 12 CPUs and 64 GiB RAM. Native tests used
the analysis environment's `bin` first on PATH, `PYTHONDONTWRITEBYTECODE=1`,
Python `-B` and pytest `-p no:cacheprovider`.

| Verification | Observed result | Scope and timing |
|---|---|---|
| Broad final workspace run | 782 passed, 78 existing warnings, 277.14 s | Exact 25 test owners below; after alias/synthetic-fixture fixes, before the final path-only `FIXTURE` export restoration |
| Exact clean candidate | 317 passed, 108.58 s | Git HEAD plus the 21 owned files; nine affected test owners; includes the restored export and no unpublished fixture data |
| Unrestricted pytest | 37 existing collection errors, 9.89 s; exit 2 | Exact same 37 error paths as baseline, with no introduced or removed error |
| Static checks | Passed | 15 Python owners compile in memory; three reviewed numerical-function ASTs unchanged; six AF defaults unchanged; no adapter edits |
| Documentation | Passed | Canonical policy synchronization; wiki build check validates 50 pages; whitespace diff check |

The broad run's warnings are existing diagnostics, not waived test failures.
The final export restoration adds only the actively imported pathname/comment;
no runtime source or test function body changed afterward. Root verified byte
identity between current source/YAML and all nine affected test owners and the
317-case tested candidate. The export does not read the local manifest; the
published batch tests generate synthetic GWAS/sample-sheet data. The separate
local compressed tests retain their existing optional fixture-path API.

Exact recorded broad invocation (`postfix-test-command.txt`), from repository
root; its listed paths are the full 25-owner run:

```bash
PATH=/Users/JJOHN41/miniconda3/envs/postgwas/bin:$PATH PYTHONPATH=src PYTHONDONTWRITEBYTECODE=1 /Users/JJOHN41/miniconda3/envs/postgwas/bin/python -B -m pytest -q -p no:cacheprovider tests/test_delimiter_header_selection.py tests/test_harmonisation_af_acceptance.py tests/test_harmonisation_build_fallback.py tests/test_harmonisation_config.py tests/test_harmonisation_exclude_chr_x.py tests/test_harmonisation_exclude_x_accounting.py tests/test_harmonisation_exclude_x_interface.py tests/test_harmonisation_exclude_x_references.py tests/test_harmonisation_exponentiated_se.py tests/test_harmonisation_exponentiated_se_integration.py tests/test_harmonisation_frequency_reference_reuse.py tests/test_harmonisation_metadata_preamble.py tests/test_harmonisation_multiallelic_reference_frequency.py tests/test_harmonisation_reference_frequency_export.py tests/test_harmonisation_reference_frequency_provenance.py tests/test_harmonisation_sample_sheet_generator.py tests/test_harmonisation_shared_coordinate_identifier.py tests/test_harmonisation_study_report_integration.py tests/test_harmonisation_replacement_af_validation.py tests/test_harmonisation_replacement_af_integration.py tests/test_harmonisation_replacement_af_reporting.py tests/test_harmonisation_eaf_replacement_resources.py tests/test_harmonisation_eaf_replacement_science.py tests/test_harmonisation_batch_preflight.py tests/test_harmonisation_resource_preflight.py --basetemp /private/tmp/postgwas-replacement-release-7390gtl3/fixtures-final-alias-fix
```

Reproducible equivalent of the clean-candidate command, from the candidate root
with candidate-only `src` on PYTHONPATH; these are all nine affected owners. The
actual run used an owned temporary `--basetemp`; its suffix was not captured in
the command artifact, so none is invented here. A fresh owned `--basetemp` may be
added when reproducing this command:

```bash
PATH=/Users/JJOHN41/miniconda3/envs/postgwas/bin:$PATH PYTHONPATH=src PYTHONDONTWRITEBYTECODE=1 /Users/JJOHN41/miniconda3/envs/postgwas/bin/python -B -m pytest -q -p no:cacheprovider tests/test_harmonisation_af_acceptance.py tests/test_harmonisation_frequency_reference_reuse.py tests/test_harmonisation_replacement_af_validation.py tests/test_harmonisation_replacement_af_integration.py tests/test_harmonisation_replacement_af_reporting.py tests/test_harmonisation_eaf_replacement_resources.py tests/test_harmonisation_eaf_replacement_science.py tests/test_harmonisation_batch_preflight.py tests/test_harmonisation_resource_preflight.py
```

Equivalent unrestricted verification, from repository root:

```bash
PATH=/Users/JJOHN41/miniconda3/envs/postgwas/bin:$PATH PYTHONPATH=src PYTHONDONTWRITEBYTECODE=1 /Users/JJOHN41/miniconda3/envs/postgwas/bin/python -B -m pytest -q -p no:cacheprovider
```

Audit evidence: `postfix-tests.log`, `candidate-verified-tests.log`,
`verified-full-pytest.log`, `postfix-static-check.txt`, and
`final-rereview.md` (including its approved export addendum) in the invocation's
`.superpowers/sdd/2026-11-04-replacement-af-validation/` directory. This permanent
record preserves their material outcomes without requiring those local artifacts
to be committed.

### Six ledger rulings and their costs

1. Reuse the active feature checkout and snapshots rather than create a new
   worktree. This preserves the user's checkout; incorrect attribution would
   require removing an unrelated staged hunk before commit.
2. Investigate pre-existing runtime failures without an unnecessary approval
   pause. Authorized in-scope fixes proceed, while unrelated failures remain
   explicit verification limits; this does not authorize unrelated repairs.
3. Reuse the existing memory budget/safety controls for exact global Spearman.
   A null direct-library budget remains unbounded. No new setting is needed,
   but an underestimated requirement can exceed available memory; estimates
   are not OS limits and no disk-spill rank path is provided.
4. Review task-owned patches against the recorded integration base and create
   one final authorized local commit rather than intermediate feature commits.
   This avoids publishing a provisional dependent contract; incorrect patch
   attribution still requires correcting staged hunks before committing.
5. Bound preparation structurally to one projected chromosome/reference, without
   treating the later whole-worker scheduler as a preparation estimator. This
   avoids a misleading estimate containing unrelated snapshots/VCF buffers and
   worker allocations; a large chromosome/reference can still exceed memory.
   Preparation admission planning requires a separate scheduling change.
6. Include focused feature/regression tests and their resource-preflight helper
   despite the broad local-test ignore. `.gitignore` remains unchanged; the cost
   is additional test source in Git, instead of dependence on unpublished local
   fixtures. No private GWAS payload or local manifest is published.

### Final scientific and provenance corrections

The approval contract retains lexical `source_path`, `staged_path` and
`reference_path` separately from resolved identities. Activation preserves
formatted declared paths before resolving aliases. The same pathname is checked
before diagnostics, after measurement and at adoption; the verified physical
attachment path is pinned before intervening strand work. Alias retargeting
cannot substitute unmeasured frequencies. This is path pinning, not a lock
against arbitrary concurrent writes to the physical target. Post-measurement
invalidation preserves metrics and issues no approval contract. Unchanged
symlink/hardlink reuse, legitimate row subsets and original invalid-A provenance
remain supported. Batch fixtures now generate valid synthetic data; their restored
path-only export serves an existing local consumer without making candidate tests
read that manifest.

### Benchmarks and interpretation limits

Synthetic numeric-only Pearson 100M/200M pair runs took 4.5046/7.5457 s;
20-ms sampled process RSS was 207,536,128/209,420,288 bytes. Separate shared merges
with 100K/1M study and 110K/1.1M reference rows took 0.2182/1.9767 s;
sampled RSS was 323,158,016/1,577,566,208 bytes, with AF matching the oracle at
1e-15. Numeric tests omit source I/O/merging; merge tests omit full harmonisation.
Production 100–200M reference/staging/full-workflow capacity, Linux execution,
true OS peak RSS, I/O throughput and spill-space needs remain unverified. No
production-scale workflow was launched. No cutoff calibration or proof of
per-row EAF semantics is claimed. Same-source evidence is reused; different
file paths do not establish independent scientific provenance.

### Explicit unresolved baseline checks

The separate `tools/docs/validate_wiki_cli.py --commands-only` check remains red
because its pre-existing extractor treats the value in `--name postgwas` as a
command start and extracts `postgwas --jobs 2` from valid installer examples in
`README.md` and `docs/wiki/getting-started/installation.md`. The installer supports
`--jobs`; neither examples nor validator were changed for this feature.

The invocation-start tracked baseline passed 578 tests. The final 782-case
workspace run includes a pre-existing tracked exponentiated-SE test importing
an ignored concordance helper; therefore it is not a fresh-clone guarantee for
the entire tracked suite. The 317-case clean candidate specifically verifies
the nine affected test owners and their published dependencies.

Unrestricted collection remains blocked by 36 archived-copy/module-name
collisions plus the existing ignored coordinates test's missing helper import.
The exact final error set equals the baseline set:

```text
ERROR output/installation_validation_20261030/source_candidate_02/tools/resource_preparation/dbsnp_alfa_latest/test_prepare_population_frequencies.py
ERROR output/installation_validation_20261030/source_candidate_02/tools/resource_preparation/dbsnp_build157_grch37/test_prepare.py
ERROR output/installation_validation_20261030/source_candidate_02/tools/resource_preparation/finngen_r13_annotations/test_prepare_info_scores.py
ERROR output/installation_validation_20261030/source_candidate_02/tools/resource_preparation/thousand_genomes_phase3/test_prepare.py
ERROR output/installation_validation_20261030/source_candidate_03/tools/resource_preparation/dbsnp_alfa_latest/test_prepare_population_frequencies.py
ERROR output/installation_validation_20261030/source_candidate_03/tools/resource_preparation/dbsnp_build157_grch37/test_prepare.py
ERROR output/installation_validation_20261030/source_candidate_03/tools/resource_preparation/finngen_r13_annotations/test_prepare_info_scores.py
ERROR output/installation_validation_20261030/source_candidate_03/tools/resource_preparation/thousand_genomes_phase3/test_prepare.py
ERROR output/installation_validation_20261030/source_candidate_04/tools/resource_preparation/dbsnp_alfa_latest/test_prepare_population_frequencies.py
ERROR output/installation_validation_20261030/source_candidate_04/tools/resource_preparation/dbsnp_build157_grch37/test_prepare.py
ERROR output/installation_validation_20261030/source_candidate_04/tools/resource_preparation/finngen_r13_annotations/test_prepare_info_scores.py
ERROR output/installation_validation_20261030/source_candidate_04/tools/resource_preparation/thousand_genomes_phase3/test_prepare.py
ERROR output/installation_validation_20261030/source_candidate_05/tools/resource_preparation/dbsnp_alfa_latest/test_prepare_population_frequencies.py
ERROR output/installation_validation_20261030/source_candidate_05/tools/resource_preparation/dbsnp_build157_grch37/test_prepare.py
ERROR output/installation_validation_20261030/source_candidate_05/tools/resource_preparation/finngen_r13_annotations/test_prepare_info_scores.py
ERROR output/installation_validation_20261030/source_candidate_05/tools/resource_preparation/thousand_genomes_phase3/test_prepare.py
ERROR output/installation_validation_20261030/source_candidate_06/tools/resource_preparation/dbsnp_alfa_latest/test_prepare_population_frequencies.py
ERROR output/installation_validation_20261030/source_candidate_06/tools/resource_preparation/dbsnp_build157_grch37/test_prepare.py
ERROR output/installation_validation_20261030/source_candidate_06/tools/resource_preparation/finngen_r13_annotations/test_prepare_info_scores.py
ERROR output/installation_validation_20261030/source_candidate_06/tools/resource_preparation/thousand_genomes_phase3/test_prepare.py
ERROR output/installation_validation_20261030/source_candidate_07/tools/resource_preparation/dbsnp_alfa_latest/test_prepare_population_frequencies.py
ERROR output/installation_validation_20261030/source_candidate_07/tools/resource_preparation/dbsnp_build157_grch37/test_prepare.py
ERROR output/installation_validation_20261030/source_candidate_07/tools/resource_preparation/finngen_r13_annotations/test_prepare_info_scores.py
ERROR output/installation_validation_20261030/source_candidate_07/tools/resource_preparation/thousand_genomes_phase3/test_prepare.py
ERROR output/installation_validation_20261030/source_initial/tools/resource_preparation/dbsnp_alfa_latest/test_prepare_population_frequencies.py
ERROR output/installation_validation_20261030/source_initial/tools/resource_preparation/dbsnp_build157_grch37/test_prepare.py
ERROR output/installation_validation_20261030/source_initial/tools/resource_preparation/finngen_r13_annotations/test_prepare_info_scores.py
ERROR output/installation_validation_20261030/source_initial/tools/resource_preparation/thousand_genomes_phase3/test_prepare.py
ERROR output/sample_sheet_metadata_validation_20261030/source/tools/resource_preparation/dbsnp_alfa_latest/test_prepare_population_frequencies.py
ERROR output/sample_sheet_metadata_validation_20261030/source/tools/resource_preparation/dbsnp_build157_grch37/test_prepare.py
ERROR output/sample_sheet_metadata_validation_20261030/source/tools/resource_preparation/finngen_r13_annotations/test_prepare_info_scores.py
ERROR output/sample_sheet_metadata_validation_20261030/source/tools/resource_preparation/thousand_genomes_phase3/test_prepare.py
ERROR tests/test_harmonisation_coordinates.py
ERROR tools/resource_preparation/dbsnp_alfa_latest/test_prepare_population_frequencies.py
ERROR tools/resource_preparation/dbsnp_build157_grch37/test_prepare.py
ERROR tools/resource_preparation/finngen_r13_annotations/test_prepare_info_scores.py
ERROR tools/resource_preparation/thousand_genomes_phase3/test_prepare.py
```

### Final compliance and local-commit record

Implementation/re-review and scoped tests are complete. Cumulative code,
configuration, test and documentation diffs were reviewed; the shared fingerprint
helper removed the identified duplication, and the final alias/data-portability
findings are resolved. No scientific defaults, numerical helper behavior or
protected adapters changed during final corrections. Scientific validation is
supported by the recorded numerical, integration, native-output and reporting
checks; remaining capacity and baseline-suite limits are explicit above.
Task-owned temporary validation files are cleaned by the owning root before
handoff. Only exact-path staging/local commit remains; no push is implied.
