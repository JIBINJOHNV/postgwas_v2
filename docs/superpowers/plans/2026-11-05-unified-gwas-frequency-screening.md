# Unified GWAS Frequency Screening Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Apply one GWAS-row effect-frequency screen to native, initially external and replacement inputs; use the requested folded-effect diagnostic only when needed; stop on a MAF-like replacement.

**Architecture:** Extend the existing frequency screen, shared allele-frequency preparation, reference matching and service orchestration. Preserve a row-bound initial A and candidate replacement B in original GWAS effect-allele direction. Reuse the existing staged references, joins, pooled statistics and approval contracts; do not create a parallel mapper or global multiallelic rewrite.

**Tech Stack:** Existing Python >=3.10 package, Polars and pytest. Use the installed PostGWAS environment for local verification; its verified Polars version is 0.20.31. Do not introduce dependencies or change supported version floors.

**Spec:** ../specs/2026-11-05-unified-gwas-frequency-screening-design.md

## Global Constraints

- Read the owning repository AGENTS.md and current cumulative diff before each change cycle. Preserve unrelated changes and untracked files.
- Keep source selection, existing configuration values, normal row QC and public interfaces unchanged except for the explicit transitions in the spec.
- Screen one effect-aligned value per original GWAS row; native/external/replacement must use the same original-EA direction and configured cutoff.
- Never overwrite the selected A, silently reinterpret external provenance, or adopt a failed replacement.
- Use min(f,1-f) for the requested diagnostic; remove multiallelic-only diagnostic exclusion without weakening strand/output safeguards or modifying complete-site attachment.
- Preserve fixed-N 0.50 diagnostic coverage and independent per-chromosome 0.80 attachment coverage defaults. Do not apply the earlier 0.45 sensitivity trial.
- Preserve source/policy/row-bound cache identities and bounded chromosome processing. Do not collect a full genome to screen it or rerun a heavy source join solely for screening.
- Do not modify src/postgwas/modules/harmonisation/adapters/.
- Implementation is in the isolated codex/unified-af-screening checkout; the original main checkout and server installation remain unchanged.

## Review Focus

- A huge panel with frequencies unlike its GWAS overlap: only merged GWAS rows may control the screen (Task 1).
- GWAS EA opposite REF/ALT, reverse strand or unresolved: native/external screens must use the same safe original-EA convention (Tasks 1-2).
- Initial external MAF support discovered after staging: replacement must acquire new resource proofs and must not adopt the initial source's cache (Task 3).
- Replacement that agrees perfectly with default but still looks MAF-like: stop before further comparison/adoption (Task 3).
- A global pass followed by a chromosome coverage/identity failure: stop under existing policy and record no final approval (Tasks 3-4).

## File map

Paths below are relative to /Users/JJOHN41/Documents/developing_software/postgwas_v2.

| File | Responsibility |
|---|---|
| src/postgwas/modules/harmonisation/shared/statistics.py | Existing frequency_maf_screen and reusable count pooling |
| src/postgwas/modules/harmonisation/study_properties.py | Native selected-frequency evidence and common screen decision |
| src/postgwas/modules/harmonisation/strand.py | Prepared original-EA reference alignment and scoped diagnostic eligibility |
| src/postgwas/modules/harmonisation/allele_frequency.py | Selected A/B, merged screening, folded diagnostics, source transitions and approval |
| src/postgwas/modules/harmonisation/service.py | Staging/activation order and hand-off to workers |
| src/postgwas/modules/harmonisation/study_reporting.py | Screen scope, evidence, stop reason and source adoption reporting |
| src/postgwas/modules/harmonisation/html_report.py; vcf_provenance.py; concordance/service.py | Existing report and downstream source/provenance consumers |
| src/postgwas/config/defaults/modules/harmonisation.yaml | Updated help and four additive VCF provenance-header mappings; preserve scientific defaults |
| docs/modules/harmonisation/eaf_replacement.md; configuration.md | User-facing flow, formulas, retained thresholds and failure behavior |
| docs/wiki/harmonisation/policies.md; processing-order.md; sample-sheet.md; outputs-and-qc.md | Public policy and processing documentation; regenerate generated content through existing tooling |
| Existing harmonisation test files | Regression and integration coverage listed below |

### Task 1: Common selected-frequency screening

**Files:** shared/statistics.py, study_properties.py, allele_frequency.py; tests/test_harmonisation_external_af_preflight.py and existing frequency-screen/study-property tests.

**Interfaces:** Preserve frequency_maf_screen(frame, column, cutoff) behavior for existing DataFrame callers. Reuse its count definition for pooled chromosome streams rather than duplicating the cutoff rule. Preserve existing public validation function signatures. Extend existing study_decisions provenance and row-bound work data to carry A without pretending an external column was native.

- [x] Write failing integration tests for merged GWAS scope, unrelated panel rows and original-EA alignment. An extra million reference-only rows must leave the screen unchanged. Unequal chromosome sizes must use pooled counts, not averaged percentages.
- [x] Pin the boundary assertions: 190/200 usable values <=0.5 at cutoff .95 is not MAF-like; 191/200 is MAF-like; AF=.5 counts toward the numerator; AF=0 and 1 are usable; null, nonnumeric, NaN, infinity and out-of-range values are unusable; zero usable values stop.
- [x] Verify the new behavior tests fail under current code for the intended reason. Run existing helper tests as boundary controls.
- [x] Implement the all-source-row, original-EA-aligned screen stream. Do not screen the later diagnostic-only working frame, which has already lost unmatched/palindromic/ambiguous rows. Unresolved safe assignments remain null with explicit counts; they must not be guessed.
- [x] Preserve A and its mapping/source identities. Keep the initial input payload unchanged. Source-only raw table statistics may remain descriptive but must not trigger this screen.
- [x] Run the focused screen, study-property and initial-external tests; inspect observed fractions and exact row identities. Review cumulative diff before committing this task through the normal development workflow.

### Task 2: Default-reference folded-effect diagnostic

**Files:** strand.py, allele_frequency.py; tests/test_harmonisation_allele_frequency.py, test_harmonisation_eaf_replacement_science.py, test_harmonisation_external_af_preflight.py and test_harmonisation_multiallelic_reference_frequency.py.

**Interfaces:** Keep prepare_study_reference_matches and the existing pooled comparison APIs. Generalize scoped diagnostic candidate selection/alignment so resolved multiallelic status alone does not reject a pair; preserve structural metadata and unchanged callers' strand/output protections.

- [x] Write failing tests showing a non-MAF-like A skips semantic default comparison, while MAF-like A enters it.
- [x] Add direct, swapped and reverse-complement cases proving the default f and A use the same original GWAS EA. Assert the requested f=.8 -> 1-f=.2 -> folded=.2 and f=.3 -> folded=.3 calculations, with unchanged A.
- [x] Add resolved multiallelic diagnostic cases, including zero-frequency extra ALT records. Assert inclusion under the requested folding convention while existing attachment outputs and final strand behavior remain unchanged.
- [x] Keep tests for unresolved orientation, palindromes, incompatible builds, missing selected AF and conflicting keys. This change must not make those silently eligible.
- [x] Implement default effect-frequency diagnostic alignment through the existing machinery. Avoid recovering a different named allele by an accidental intermediate source-ALT complement. Scope the 1-f diagnostic convention separately from complete-site final assignment.
- [x] Reuse the current dual-rule comparator, counts, fixed N, Pearson/Spearman behavior and configured thresholds. Test exact boundary values and both-rule/neither-rule failures. Retain the independent-source requirement for initially external A.
- [x] Keep raw aligned r and a separate folded view. Add an assertion that the unchanged major-fraction rule still sees raw r>0.5; globally folded reference storage must not silently make replacement selection impossible.
- [x] Update tests whose old contract required multiallelic-only diagnostic exclusion; preserve their attachment assertions. Run the focused diagnostic suites and inspect the cumulative diff before committing this task.

### Task 3: One replacement transition and strict merged-B screen

**Files:** allele_frequency.py, service.py; tests/test_harmonisation_external_af_preflight.py, test_harmonisation_replacement_af_integration.py, test_harmonisation_replacement_af_validation.py, test_harmonisation_eaf_replacement_resources.py, test_harmonisation_eaf_replacement_science.py and test_harmonisation_external_reference_staging.py.

**Interfaces:** Retain existing callers of validate_replacement_allele_frequency and merge_external_allele_frequencies; add an optional staging_summary keyword to replacement validation for producer proof. Generalize selected-A provenance so A may come from native or initially external frequency. Reuse _activate_eaf_replacement_resources, staging and existing identity/caching contracts when the transition occurs.

- [x] Write failing tests for native-A -> MAF support -> replacement and external-A -> MAF support -> replacement. Require both public options and canonical YAML paths when replacement is needed; do not select replacement for an accepted A.
- [x] Write a failing test that merges B, screens its GWAS-row original-EA values and stops when MAF-like. Assert no subsequent default comparison, recursive replacement, fallback, chromosome worker or final approval occurs, even if a default comparison would be numerically excellent.
- [x] Test no usable B, cutoff equality, unequal chromosome sizes, different raw-panel versus merged-row distributions, direct/swapped keys, missing matches and invalid values.
- [x] Implement the resource lifecycle. Initial external MAF support is discovered after initial staging, so activate/rebind/stage the selected replacement at that point. Preserve frozen A and its source proof; never reuse the original external assignment descriptor as replacement approval. Native replacement selection may retain its existing earlier activation where equivalent.
- [x] Generalize native-era state assumptions consistently in _activate_eaf_replacement_resources, selected_eaf_replacement, _require_replacement_validation and worker adoption. External A must remain private row-bound evidence with a descriptor; do not manufacture a native eaf_col. Preserve original-EA raw B separately from its folded view.
- [x] Screen B before the remaining replacement consistency check. A MAF-like result is terminal for that dataset and names the replacement file/column and observed cutoff evidence. No automatic retry/fallback is permitted.
- [x] If B passes, retain the existing A-versus-min(B,1-B) consistency thresholds and the new diagnostic eligibility convention, then publish only unchanged, approved assignments. Preserve existing same-source reuse labels and identity checks.
- [x] Update the old test_declared_replacement_does_not_screen_its_minor_alt_values to assert the newly requested behavior rather than preserving the obsolete exemption. Update old replacement diagnostic-exclusion expectations without weakening attachment/missing-value tests.
- [x] Verify source mutation, alias changes, staged-content changes, partial publication, source-name collisions and worker-subset/reordering cases. Test that a global pass cannot override per-chromosome coverage failure. Run focused integration suites and review the complete diff before committing this task.

### Task 4: Reporting, documentation and end-to-end verification

**Files:** study_properties.py source-selection/finalization, study_reporting.py, html_report.py, vcf_provenance.py, concordance/service.py, canonical YAML help text, harmonisation and wiki documentation; tests/test_harmonisation_replacement_af_reporting.py and existing report/provenance/service tests.

**Interfaces:** Extend existing evidence records consumed by terminal, canonical logs, manifests and HTML. Use the existing evidence/reporting pipeline; do not independently reanalyse frequencies while rendering reports.

- [x] Add failing reporting tests for initial selected A, original-EA screen scope, skipped default comparison, selected B, rejected MAF-like B and successful final source adoption. No unusable/failed input may be labelled accepted.
- [x] Record the required source/column/identity, row counts, screen numerator/denominator/cutoff, fixed diagnostic N/n, comparison and replacement results, exclusions and exact failure reasons. Preserve the distinction between supplied external AF and study-measured AF provenance.
- [x] Update current documentation and configuration help to describe the common GWAS-row screen, scoped 1-f convention and strict replacement rejection. Keep all numeric defaults unchanged and avoid claiming the heuristic proves MAF semantics.
- [x] Verify downstream concordance/export provenance retains absent-native-column status while identifying A's external source and the final approved B. Use existing documentation generators for generated policy tables.
- [x] Run representative suites using the installed environment, for example: `PYTHONDONTWRITEBYTECODE=1 /Users/JJOHN41/miniconda3/envs/postgwas/bin/python -B -m pytest -q -p no:cacheprovider tests/test_harmonisation_external_af_preflight.py tests/test_harmonisation_replacement_af_integration.py tests/test_harmonisation_replacement_af_validation.py tests/test_harmonisation_eaf_replacement_science.py tests/test_harmonisation_multiallelic_reference_frequency.py tests/test_harmonisation_strand.py`.
- [x] Run the repository's full required test/check commands. Report all failures, skips, unavailable dependencies and unverified behavior by name; do not claim readiness from focused tests alone.
- [ ] Replay the available real fetal-WLM extract at unchanged thresholds, comparing initial/merged screens, N/n, mapped direction, output AF and stop/adoption states. Label it representative rather than the exact server pair file.
- [ ] Validate the exact full server wlm_sex_gest_plus_p_fetal_wlm_pair.gz with its resolved YAML/resources before claiming the reported run is fixed. Capture a complete outcome through chromosome coverage and final output checks. If server access is unavailable, mark that validation outstanding.
- [x] Obtain an independent review of scientific alignment, state transitions, cache/source integrity and cumulative diff. Remove only generated task-owned temporary files. Apply the owning AGENTS.md compliance checklist and report the remaining limitations before handoff.

## Acceptance criteria

The implementation is ready for review only when native/external A share the same GWAS-row screen, default diagnostics run only for a MAF-like A, requested folded-effect calculations use the same EA direction, resolved multiallelic status alone does not exclude diagnostic evidence, replacement B is screened and rejected immediately if MAF-like, retained consistency/normal QC still applies, and all source identities/provenance remain traceable. Passing representative local tests is not a claim of full server success.

## Execution notes

The requested dataset-only decisions and cached worker adoption are implemented. Numeric defaults remain unchanged. Native multiallelic output swaps now reject/fail under the existing strand ambiguity policy because 1-native-EAF cannot recover a particular ALT AF; diagnostic folding and external complete-site attachment are separately preserved. Four additive VCF metadata mappings retain initial external A identity and absent native AF. Final stable regression status is recorded in the linked audit report. The available real chr22 check covers AF preparation only; exact server-pair/full-output validation remains outstanding.

Final local verification: 1,200 selected tests passed (37 files), policy documentation synchronized, wiki 50-page check passed. Unrestricted local collection, unchanged historical baseline comparisons and two unchanged CLI documentation extraction findings remain documented pre-existing limitations. The representative real chr22 AF-preparation replay passed; exact server/full-output validation remains unchecked. See `/Users/JJOHN41/Downloads/gwas_sumstat/harmonisation_mapping_audit_20261105/implementation_validation.md`.
