# Whole-study AF acceptance implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Require whole-study usable reference overlap of at least `0.50`, correlation `>= 0.90`, and aligned reference-minor fraction `> 0.90`; preserve the existing shared directional error margin at `0.02`.

**Architecture:** Reuse current whole-study pair preparation and extend `_compare_frequency_pairs(..., accept_direct_eaf=True)` with the study total. Add overlap/reference-minor-fraction policies in canonical YAML, reuse the existing shared correlation policy at `0.90`, preserve the shared directional margin at `0.02`, and propagate counts and thresholds through evidence and reporting.

**Tech Stack:** Python, Polars, canonical YAML policy registry, existing configuration exporter, pytest.

**Spec:** `docs/superpowers/specs/2026-11-04-whole-study-af-acceptance-design.md` (approved by the user for implementation and commit, with the error-margin change explicitly excluded).

## Global constraints

- Read and apply `AGENTS.md` before every change cycle; inspect the cumulative diff and preserve overlapping user work.
- Correlation `>= 0.90`; reference fraction `> 0.90`; numerator uses the existing `reference AF <= 0.5` definition.
- Require at least `0.50` usable comparisons relative to all rows entering whole-file AF preflight, retaining the existing `1,000` absolute floor and `0.10` mean-difference maximum. The new fraction gates apply only to the current whole-study MAF-like native-AF path.
- One canonical schema-valued YAML source for thresholds/defaults; no Python defaults, duplicate configuration files or new matching implementation.
- Change existing shared `eaf.maf_reference_min_correlation` to `0.90`, including its existing external/folded-MAF uses. Preserve `eaf.maf_reference_error_margin: 0.02`. Preserve source selection, original AF values, chromosome QC, absolute overlap semantics and protected adapters.
- No additional study collection, new scientific calibration claim, or 100–200-million-row performance claim.

## Review focus

- The shared correlation intentionally tightens ordinary external/folded-MAF paths; the error margin remains unchanged; no redundant direct-only correlation key should remain in the proposed interface.
- Exactly `0.90` must pass correlation but fail the reference fraction gate; exactly `0.5` remains included in the fraction numerator.
- The overlap denominator must remain all study rows entering preflight, not matched or eligible rows; test fractional equality, odd counts and the absolute floor without a floating-point ceiling error.
- Failing only the new fraction gate must not imply MAF, select replacement, or allow weaker fallback acceptance.
- Undefined constant-vector correlations must be reported as unavailable even if the installed backend returns a finite artifact; preserve valid nonconstant-vector decisions.
- Exported settings, recorded direct thresholds, checkpoint identity and visible failure explanations must agree.
- These gates can accept a true-MAF control when most effect alleles are minor; retain an explicit methodological limitation.

## Task 1: Implement the configured whole-study acceptance rule

**Files:**

- Modify `src/postgwas/config/defaults/modules/harmonisation.yaml`.
- Modify `src/postgwas/config/models/modules/harmonisation.py` to reuse registry validation at the shared configuration boundary while preserving raw overrides.
- Modify `src/postgwas/modules/harmonisation/allele_frequency.py`.
- Modify `src/postgwas/modules/harmonisation/study_reporting.py`.
- Modify `docs/modules/harmonisation/configuration.md`, `docs/modules/harmonisation/eaf_replacement.md`, and `docs/wiki/harmonisation/policies.md` where they describe this decision.
- Add self-contained committed regressions in `tests/test_harmonisation_af_acceptance.py`. Verify the existing local `tests/test_harmonisation_study_eaf.py`, `tests/test_harmonisation_config_export.py`, `tests/test_harmonisation_af_failure_reporting.py`, `tests/test_harmonisation_study_reporting.py`, and `tests/test_harmonisation_eaf_replacement_science.py` without publishing their ignored historical suite or real-data fixtures.

**Interfaces:**

- Consumes existing `Policies.get`, valid pooled pair columns, `accept_direct_eaf`, and the existing dataset decision structure.
- Produces new `eaf.maf_reference_min_overlap_fraction: 0.50` and `eaf.study_reference_min_minor_fraction: 0.90`; changes existing `eaf.maf_reference_min_correlation: 0.90` and preserves `eaf.maf_reference_error_margin: 0.02`. Retires the never-implemented direct-only correlation proposal.
- Extends `_compare_frequency_pairs` with optional keyword `study_variant_count=None`, retaining the return type and existing external-call behavior apart from explicitly changed shared correlation default. Whole-study calls require a positive total supplied once from `df.height`.
- Adds evidence fields `study_variant_count`, `observed_overlap_fraction`, `minimum_overlap_fraction`, and `minimum_reference_effect_allele_minor_fraction`; keeps existing `minimum_overlap` as the absolute floor.
- Preserves `reference_effect_allele_minor_fraction` and uses `reference_minor_fraction_not_above_minimum` when otherwise passing direct agreement fails the new strict fraction gate.

- [x] Add failing configuration tests for new overlap `0.50`, new reference-minor fraction `0.90`, changed shared correlation `0.90`, unchanged error margin `0.02`, retained absolute count `1,000`, and retained mean-difference maximum `0.10`. Cover overrides and export round-trips in all supported styles. Reject negative, above-one, non-finite, boolean and `90` values for proportions; retain integer validation for the old count option.
- [x] Add scientific boundary tests: `test_direct_correlation_uses_shared_threshold`, `test_direct_correlation_threshold_is_inclusive`, `test_reference_minor_fraction_requires_strictly_more_than_threshold`, and `test_reference_half_frequency_counts_in_minor_fraction`. Use finite paired vectors, exact reference-minor count ratios `900/1000` and `901/1000`, and measured-correlation equality through a threshold override to test the inclusive operator without an undocumented numerical tolerance.
- [x] Add `test_whole_study_overlap_uses_all_input_rows`, `test_overlap_fraction_equality_and_odd_totals`, `test_overlap_fraction_override_avoids_ceiling_roundoff`, and `test_absolute_overlap_floor_is_retained`. Test `500/1000` versus `499/1000` with a lowered absolute floor to isolate the fraction, `14/100` under a `0.14` override, and totals below `2,000` where the retained `1,000` floor dominates. Assert excluded/missing/palindromic rows remain in the denominator, the same denominator is used across all chromosomes, and insufficient coverage cannot trigger MAF replacement.
- [x] Add `test_fraction_gate_alone_cannot_select_replacement`: otherwise passing direct agreement with a failing reference fraction returns inconclusive and preserves input values, even when a replacement is configured. Assert a clear failed-gate reason. Retain separate tests for independently supported MAF and missing replacement arguments.
- [x] Preserve the existing constant-AF failure-reporting expectation and add direct/folded constant-vector cases for supplied and reference values, under both Pearson and Spearman. Confirm unavailable correlations are represented as `None` rather than finite artifacts.
- [x] Run the five Task 1 test files and confirm new tests fail because the requested policies/gate do not exist yet. Record baseline failures separately from regressions.
- [x] Add the two new canonical YAML proportion policies and update the existing shared correlation default while preserving the error-margin default. Include validator metadata, scope, operators, unchanged-margin explanation, scientific caveat, EAF key ordering and common export. Reuse the registry loader, preserve integer overlap semantics, and do not add duplicate Python defaults or the superseded direct-only correlation key.
- [x] Verify the actual public `config validate` command for full-run and module-only exports. Reuse `load_policies` at the shared model boundary, translate `PolicyError` to the existing configuration-error path, reject both invalid overrides together, preserve nested/dotted numeric-string overrides and export/checkpoint identity, and retain base configuration-only operation without Polars imports.
- [x] Update `allele_frequency.POLICY_KEYS`. Extend `_compare_frequency_pairs` in place with `study_variant_count=None`; require a positive total on the whole-study path, calculate `paired.height / study_variant_count`, enforce both count and fraction gates, and reuse the existing reference-minor fraction metric with its strict minimum. A failed overlap fraction or a failed minor-fraction gate alone returns inconclusive and cannot authorize replacement. Existing external calls omit the total and retain count-based overlap while inheriting the explicitly changed correlation and unchanged margin.
- [x] Guard direct and folded correlation expressions when either paired vector has equal minimum and maximum, following `src/postgwas/modules/harmonisation/concordance/analysis.py:889`. Keep the check in the existing comparator, use no numerical epsilon or uniqueness set, and preserve valid nonconstant-vector results.
- [x] Update `validate_study_allele_frequency`: pass the single whole-study `df.height` denominator, preserve pooled matches and source rows, and propagate absolute overlap, observed/minimum fractions, shared correlation threshold and actionable failure evidence. Continue the separate MAF-support path only when coverage is sufficient and direct agreement itself fails.
- [x] Update `study_reporting.py` and existing failure reporting to show observed versus required values and exact operators. Tests must cover pass, fraction failure, correlation failure and independent MAF replacement, with no misleading `.80` direct threshold.
- [x] Update method/configuration documentation with the full run-YAML keys, inclusive `AF <= 0.5` definition, strict reference-minor fraction, inclusive overlap, fixed denominator, retained absolute floor and mean-error gate, shared threshold scope and synthetic-control limitation. Explain that the existing `0.02` error margin is unchanged and is not a per-variant tolerance or maximum error.
- [x] Run the Task 1 tests, inspect the complete cumulative diff and confirm all pass without weakening scientifically valid tests. Update expected shared threshold values while preserving independent coverage of direct versus folded correlation calculations and the distinction between the new lower reference-minor bound and existing upper directional bound.

## Task 2: Verify downstream behavior and configuration identity

**Files:**

- Extend the self-contained `tests/test_harmonisation_af_acceptance.py` where integration coverage is needed; verify existing local `tests/test_harmonisation_allele_frequency.py`, `tests/test_harmonisation_eaf_replacement_downstream.py`, `tests/test_harmonisation_eaf_replacement_interface.py`, and `tests/test_harmonisation_policy_defaults.py`.
- Update the tracked `tests/test_harmonisation_config.py` expectation for the deliberately changed shared correlation default; retain the separate non-effect-frequency default of `0.80`.
- Extend existing applicable tests only where needed; production changes beyond Task 1 require a demonstrated integration dependency and cumulative diff review.

**Interfaces:**

- Consumes Task 1 policy/evidence fields and the unchanged prepared-match/worker interface.
- Produces regression evidence that the new policy is enforced once per study and participates in existing checkpoint identity without changing other source paths.

- [x] Add/update regression cases proving accepted minor-effect-allele EAF remains unchanged during preflight, prepared reference matches are reused, chromosome allele swaps still occur once, and existing replacement provenance is preserved.
- [x] Verify ordinary external/folded-MAF cases use the changed shared correlation `0.90` and unchanged margin `0.02`. Correlation exactly `0.90` passes; values between `0.80` and `0.90` no longer pass the default correlation gate. Error gap exactly `0.02` is inconclusive and a strictly larger gap passes when all other gates pass. Directly accepted whole-study EAF remains independent of the directional error-margin gate.
- [x] Add a deterministic true-MAF control with 9,400 minor-reference values spanning `0.01..0.49` and 600 major-reference values spanning `0.55..0.65`; supplied AF is `min(reference, 1-reference)`. Verify the measured correlation is approximately `0.9523683946`, minor fraction `0.94`, and mean error `0.012`. Record that these user-selected gates pass this control; do not claim semantic identification or invent an extra rejection rule.
- [x] Verify no-overlap, constant/non-finite vectors, invalid source values, excluded diagnostic sites, the global denominator and chromosome pooling remain covered. Keep external attachment's separate `eaf.external_min_match_fraction: 0.80` unchanged; preserve stop behavior for unresolved decisions.
- [x] Verify exported or overridden policy changes alter the resolved scientific configuration/checkpoint identity, following existing resume/restart policy. Reuse the shared checkpoint tests rather than adding a new checkpoint implementation.
- [x] Run the focused verification command below, review the full cumulative diff, remove only task-created caches/artifacts and provide the required change-cycle compliance report. No production-ready or full-scale claim is permitted without the corresponding evidence.

## Focused verification command

Run in the existing PostGWAS environment, with bytecode/cache disabled and `--basetemp` pointing to a disposable temporary directory:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python -B -m pytest -q -p no:cacheprovider \
  tests/test_harmonisation_af_acceptance.py \
  tests/test_harmonisation_study_eaf.py \
  tests/test_harmonisation_config_export.py \
  tests/test_harmonisation_af_failure_reporting.py \
  tests/test_harmonisation_study_reporting.py \
  tests/test_harmonisation_eaf_replacement_science.py \
  tests/test_harmonisation_allele_frequency.py \
  tests/test_harmonisation_eaf_replacement_downstream.py \
  tests/test_harmonisation_eaf_replacement_interface.py \
  tests/test_harmonisation_policy_defaults.py
```

Expected: all focused tests pass, whole-file overlap and minor-fraction gates use the specified denominators/operators, direct and folded comparisons use the requested shared correlation, directional classification preserves the existing margin, and configured values agree with exports and reporting. Run any applicable repository checks discovered during implementation; report checks that cannot run and their remaining risk.

## Planning review record

The current code, canonical registry, export loader, reports and representative tests were inspected. This revision records the user's shared correlation change and explicit exclusion of the error-margin change and the new whole-file usable-overlap fraction, preserving the inclusive half-frequency definition and the distinction between minimum reference-minor acceptance and the existing maximum for directional evidence. The source-specific whole-file/chromosome asymmetry, unchanged WLM error-margin failure and true-MAF counterexample remain explicit limitations. This plan does not authorize silently moving external attachment or refactoring full-scale input processing.

Initial baseline verification on the old implementation: `tests/test_harmonisation_config_export.py`, `tests/test_harmonisation_study_reporting.py` and `tests/test_harmonisation_af_failure_reporting.py` produced `54 passed, 1 failed`. The implemented constant-vector guard corrects that diagnostic without weakening the existing scientific test. A later real public-CLI check exposed invalid policies being accepted by configuration-only validation; the shared registry-backed model validator now rejects them through the existing `ConfigurationError` boundary. One tracked exclusion test now expects that earlier exact error class; its six invalid inputs and diagnostic checks remain unchanged.

The committed feature regressions are self-contained and do not depend on the ignored historical suite or its real-data fixtures. Optional local historical orientation/preparation cases use explicit coverage settings to isolate their original behavior, while the new default coverage rejection is tested separately. Full production GWAS execution, empirical threshold calibration and 100–200-million-row performance validation remain unperformed.


## Final verification record — 2026-11-04

- Implemented both plan tasks and the necessary shared public-configuration validation correction. The error margin remains `0.02`.
- The committed standalone AF file contains 84 self-contained cases. The final ten-file AF run (the command above, in the existing PostGWAS environment with temporary basetemp) passed **292 tests** in 141.84 seconds; 81 existing Polars coalesce deprecation warnings came from unchanged concordance code.
- A final explicit run of all 14 Git-tracked test files passed **453 tests** in 203.03 seconds; 42 existing concordance deprecation warnings remain. The focused and tracked sets overlap; these counts are not additive.
- The shared configuration/registry/export correction passed **221 tests**. Actual public CLI checks reject both invalid policy values together with exit `2`, without a traceback; valid exports, raw override shape and base operation without Polars imports are covered.
- Bare `pytest` was attempted on the final production source and stopped at the same **37 pre-existing collection errors**: copied output/tool test module collisions and `tests/test_harmonisation_coordinates.py` importing the unavailable `_stage_indexed_csv_batches`. The unrestricted suite could not execute. Those unrelated files were not changed.
- Independent Task 1 and final whole-change reviews passed; the shared-validator fix and its exact-exception test adaptation also passed scoped review. Python syntax, canonical YAML defaults, local document references, protected-path checks and cumulative whitespace checks passed.
- Scope checklist: passed for the implemented change; cumulative diff reviewed; hardcoding/duplication audit passed; controlled scientific validation completed. Global collection blockers, threshold-calibration limits and unperformed production/100–200-million-row runs remain explicitly outstanding.
