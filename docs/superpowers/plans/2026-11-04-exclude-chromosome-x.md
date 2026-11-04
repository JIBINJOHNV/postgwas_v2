# Exclude Chromosome X Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox syntax for tracking.

**Goal:** Add one explicit CLI flag that excludes chromosome X throughout harmonisation, while retaining X under the packaged defaults.

**Architecture:** Resolve the flag into the existing chromosome policy registry. Reuse chromosome-label normalization and rejection accounting to remove canonical X rows before scientific QC or study-wide decisions, and use the same exclusion policy for resource planning and validation.

**Tech Stack:** Python, argparse, canonical YAML/policy validation, Polars, existing PostGWAS logging/checkpoint infrastructure, pytest.

**Spec:** The user request and clarification are captured in the design contract below. Implementation and commit were authorized on 2026-11-04; execution results are recorded below.

## Design contract

- Public option: `--exclude-chr-x`, a boolean flag taking no value.
- Canonical key: `modules.harmonisation.policies.chromosome.exclude_x`; packaged default `false`.
- With the packaged configuration and no flag, existing chromosome-X handling is preserved. Explicit YAML `true` also requests exclusion; an omitted CLI flag must not erase YAML.
- When exclusion is enabled, identify X using the resolved chromosome normalization and rename map. Under the packaged map this includes `X`, `chrX`, `23`, `25`, `XY`, `PAR1`, and `PAR2`.
- Remove X immediately after chromosome identity can be resolved from separate or combined coordinate fields and after immutable source-row identities exist. X must not contribute to missing-value QC, allele/statistic validation, duplicate decisions, genome-build inference, strand consensus, study-property/MAF decisions, EAF comparison, reconstruction, partitions, workers, scientific QC denominators, or output variants.
- Reading a mixed-chromosome source and counting excluded rows remain necessary. Preserve original row provenance and raw-input accounting; do not edit the source file.
- Do not require or stage X-specific reference files, including when the configured startup reference inventory explicitly includes X. Shared whole-genome resources can still be required for retained chromosomes.
- Record exclusion as an explicit user-requested chromosome exclusion, distinct from invalid or unsupported chromosomes. Reconcile raw rows, exclusions, other rejections, and retained rows exactly once.
- If no variants remain, stop that dataset with an actionable message and the existing failure/reporting lifecycle; do not emit a successful empty scientific result.
- Changing this scientific scope must invalidate incompatible completed checkpoints under the global resume policy. A subsequent identical exclusion run must be reusable.
- Scope is the existing direct harmonisation command, including its optional validation. Harmonisation is currently registered with `pipeline_enabled=False`; this task does not add a pipeline adapter. Later analyses using its variant outputs will receive the retained scope.

## Global constraints

- Follow `AGENTS.md` for every change cycle and preserve unrelated `_to_delete/` and `tmp/` content.
- Do not edit `src/postgwas/modules/harmonisation/adapters/`.
- Keep one canonical policy and one shared exclusion implementation; do not hardcode an autosome replacement list or independent argparse defaults.
- Preserve configured chromosome order and other allowed-chromosome choices. Explicit X exclusion takes precedence over lists that include X.
- Reuse shared help, progress, logging, rejection, configuration, and checkpoint mechanisms.
- This is an explicit analysis-scope choice; it does not infer missing EAF or change EAF thresholds, allele conventions, or scientific calculations for retained variants.

## Review focus

- X aliases and combined coordinates must resolve consistently with existing chromosome policies.
- Missing or malformed scientific values confined to excluded X rows must not affect retained-variant analysis.
- An X-only dataset must fail cleanly before build inference, resource analysis, or worker launch.
- Missing X reference files must not block an exclusion run; missing required autosomal resources must still block it.
- Changing exclusion state must not reuse results produced under the opposite state, while unchanged state must resume safely.

## Task 1: Canonical configuration and public option

**Files:** `src/postgwas/config/defaults/modules/harmonisation.yaml`; `src/postgwas/modules/harmonisation/cli.py`; existing configuration/help tests.

**Interfaces:** `get_harmonisation_parser()` exposes `--exclude-chr-x`; `CLI_OVERRIDE_PATHS` maps its destination to `modules.harmonisation.policies.chromosome.exclude_x`; the existing resolver and policy registry own validation and defaults.

- [ ] Add failing cases to `tests/test_harmonisation_cli.py`, `tests/test_harmonisation_config.py`, and `tests/test_harmonisation_config_export.py`: packaged false; absent CLI attribute; flag resolves true; YAML true survives omitted flag; invalid non-boolean YAML fails; export/reload preserves the policy.
- [ ] Run the new cases and verify they fail for the absent policy/option.
- [ ] Add the canonical boolean policy with the existing flag validator and relevant policy/export metadata. Register argparse `store_true` with a suppressed default and the shared default-help formatter.
- [ ] Extend `tests/test_cli_help.py` for the new option, canonical `Default:` label, narrow/wide alignment, and a complete parseable example.
- [ ] Run these focused files with bytecode/cache disabled; inspect the cumulative diff and complete the repository compliance gate.

## Task 2: Early exclusion, resource planning, and complete accounting

**Files:** `src/postgwas/modules/harmonisation/shared/variant_columns.py`; `coordinates.py`; `summary_statistics_io.py`; `service.py`; `rejects.py`; affected QC/reporting and optional concordance consumers within that package.

**Interfaces:** Extend the existing shared chromosome-scope calculation to honor `chromosome.exclude_x` while preserving configured order. Reuse/refactor the current coordinate extraction and normalization path for the early filter instead of implementing a second chromosome parser. Keep the public `read_summary_statistics()` return contract stable unless an unavoidable change updates every caller together.

- [ ] Add failing fixture tests in `tests/test_harmonisation_coordinates.py` and a focused `tests/test_harmonisation_exclude_chr_x.py`: mixed autosome/X input; packaged aliases; combined coordinates; custom rename map; no-X input; missing/uninterpretable chromosome labels; X-only input.
- [ ] Assert default retention, explicit exclusion, unchanged retained values/order, immutable original rows, and disjoint exclusion/rejection totals. Assert malformed scientific values only on excluded X rows cannot trigger retained-data QC failures.
- [ ] Apply the early filter before the existing read-stage missingness gate. Provide an explicit exclusion reason/count through existing rejection and reporting infrastructure, updating every dependent reconciliation/report consumer in the same cycle.
- [ ] Use the same scope in configured startup preflight, observed-chromosome resource planning, EAF/INFO staging, worker/merge planning, and optional concordance validation. Preserve structural checks for shared files used by retained chromosomes.
- [ ] Add resource tests to `tests/test_harmonisation_resource_preflight.py`: a configured X inventory with missing X resources is ignored when excluded; required retained-chromosome failures remain errors.
- [ ] Add orchestration spies asserting X never reaches build/strand/property/AF decisions, partitioning, workers, or variant outputs. Verify X-only failure logs and progress remain below successful completion.
- [ ] Run the new fixture tests and affected coordinate, resource-preflight, input-accounting, report, and concordance suites; inspect the cumulative diff and complete the compliance gate.

## Task 3: Resume, user documentation, and final verification

**Files:** focused public-boundary tests following `tests/test_harmonisation_eaf_replacement_interface.py`; `docs/modules/harmonisation/configuration.md`; `docs/wiki/harmonisation/policies.md`; `docs/wiki/harmonisation/processing-order.md`; `docs/wiki/harmonisation/outputs-and-qc.md`; `examples/configs/harmonisation/run_config.yaml`; shared help examples as required by their existing owner.

**Interfaces:** Keep `resolve_harmonisation_configuration()` as the common effective configuration source and use the existing direct scientific checkpoint digest. Do not create a parallel resume mechanism.

- [ ] Add public-boundary tests for default-to-exclusion and exclusion-to-default checkpoint invalidation, safe identical-exclusion resume, and existing externally modified-output refusal. Confirm no X-specific dependency is newly required solely by exclusion.
- [ ] Document the default, YAML equivalent, alias behavior, filtering order, exclusions in counts/reports, X-only failure, and checkpoint consequences. Keep examples explicit about required inputs and current direct-command support.
- [ ] Proposed example after implementation:

  ```bash
  postgwas harmonisation \
    --sample-sheet /path/to/studies.csv \
    --resource-directory /path/to/resources \
    --output-directory /path/to/results \
    --exclude-chr-x
  ```

- [ ] Run all focused files named above, applicable registered-help/configuration/checkpoint regressions, the repository's relevant required checks, and `git diff --check`. Use the installed PostGWAS Python environment; record exact commands and observed results.
- [ ] Compare a representative mixed-chromosome fixture under default/exclusion modes: default preserves current X handling, exclusion passes only the expected retained rows to scientific stages, and all counts reconcile. State explicitly whether a full native GWAS run was possible.
- [ ] Review every changed call path and the complete cumulative diff; audit duplication/configuration ownership; remove task-owned temporary files; report remaining limits. Commit/push only when separately authorized.

## Planning validation

Current baseline inspected: local HEAD `774c8ca`. The direct CLI/configuration resolver, early read/coordinate ordering, dataset preparation and resource preflight, rejection accounting, checkpoint configuration integration, and existing test patterns were inspected. No runtime code or configuration was changed while preparing this plan. Implementation tests remain to be written and executed; their success is not claimed here.

## Implementation record — 2026-11-04

- Tasks 1–3 are implemented. New focused regression files cover the interfaces, early filtering, reference/concordance paths, and accounting; existing broad tests were retained.
- Independent review identified custom combined-coordinate aliases containing separators and streaming reference schema inference as two important edge cases. Both were reproduced with failing fixtures and corrected, including exclusion-aware shared EAF/INFO reads and audited reference-scope counts.
- The final focused command, using the installed PostGWAS environment with native bcftools available, passed **384 tests** with no skips. The run included all four new regression files and affected allele-frequency, strand, INFO, reference staging/preflight, rejection, coordinate, preparation, file-validation, and reporting checks.
- Existing CLI/configuration/export/checkpoint regressions separately passed 153 tests; reporting/completeness/manifest checks separately passed 174 tests. These runs overlap the final suite and are not additive totals.
- Public harmonisation help passed checks at widths 80 and 120. Wiki sources (50 pages), canonical policy-document synchronization, Python syntax checks, and diff whitespace checks passed.
- Full repository pytest remains blocked by 37 existing collection errors: copied resource tests inside ignored output directories collide with actual tool-test modules, and a local ignored coordinate test imports an obsolete staging helper. Existing broader help tests also flag unchanged MAGMA wording and ignored copied README content. The documented-command checker flags unchanged installer `--jobs` examples. These unrelated files were not modified to conceal those failures.
- The four new regression files are included with the change despite the general local-tests ignore rule. No adapters, user inputs, `_to_delete/`, or `tmp/` content was changed.
- A full server-scale harmonisation run was not performed. Existing global resource-root fingerprinting remains conservative: files under that root may remain checkpoint dependencies even if their chromosome is excluded; declared replacement templates additionally filter explicit X matches.

The original checkboxes above preserve the implementation recipe. This execution record supersedes their planning-time unchecked state and distinguishes verified feature behavior from the outstanding repository-wide checks.
