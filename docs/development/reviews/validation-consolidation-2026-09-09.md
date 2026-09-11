# Combined file validation — implementation and verification

Scope: shared file-validation work and presentation in public direct commands
and pipelines. This is not a full-data, all-engine execution or installation audit.

## Behaviour

- Reuse successful, read-only checks only for the same file identities and
  validator contract within one invocation, including nested dependencies.
  Direct reuse does not activate pipeline-only requirements. Changed or missing
  direct dependencies reach the original validator again; changed pipeline
  inputs remain an error. Failures are not cached.
- Combine each exact file path into one version-2 validation-audit entry.
  Preserve distinct check observations, consumers, warnings, failures and
  conflicting metric values. Do not replace scientific checks with availability.
- Show established file metadata once. Combine explicitly file-linked module
  fields with that file's section. Later findings use its stable file number
  without repeating the filename or prior facts. Cross-file outcomes retain
  their meaning and refer to file numbers.
- Keep the same presentation in the durable screen log under `--hide-screen`.
  Combine structured file-validation log details in the audit, while preserving
  scientific step values for downstream consumers.
- Preserve native diagnostics, commands, analysis results and progress events.
  Do not suppress arbitrary repeated text or guess ambiguous file attribution.
- Reuse exact Manhattan runtime/header probes and FLAMES import probes;
  retain native probe diagnostics. Remove the enrichment dispatcher's duplicate
  input-gene count; its analysis service remains responsible for that count.

Configuration addition: `logging.file_validation.outcome_metric_aliases`, defined
in application YAML and schema-validated. It maps presentation labels/keys to
existing metric vocabulary; a target list selects recorded evidence before its
fallback. An ordinary record count must not be labelled an index count without
the corresponding evidence. Conflicting counts remain distinct. No scientific
threshold, transformation, filtering rule or CLI default was changed.

## Files changed in this task

Paths below are relative to the repository, not an inventory of unrelated
pre-existing working-tree changes.

| Area | Files |
|---|---|
| Shared evidence and presentation | `src/postgwas/core/input_validation.py`, `src/postgwas/core/validation_reporting.py`, `src/postgwas/core/pipeline_logging.py`, `src/postgwas/core/ui/progress.py` |
| Pipeline integration | `src/postgwas/pipeline/cli.py`, `src/postgwas/pipeline/validation_reporting.py` |
| Configuration | `src/postgwas/config/defaults/application.yaml`, `src/postgwas/config/models/logging.py` |
| Module integrations | `src/postgwas/modules/filtering/reporting.py`, `src/postgwas/modules/filtering/sumstat_filter.py`, `src/postgwas/modules/qc_summary/reporting.py`, `src/postgwas/modules/manhattan/service.py`, `src/postgwas/modules/flames/service.py`, `src/postgwas/modules/enrichment/main.py` |
| Regression tests | `tests/test_validation_consolidation.py` (new), `tests/test_file_validation_reporting.py`, `tests/test_pipeline_validation_reporting.py` |
| Rules and documentation | `AGENTS.md`, `docs/wiki/core/pipeline-input-validation.md`, `docs/wiki/core/pipeline-workflow.md`, this verification record |

The superseded pipeline-specific print helper was removed in favour of the
shared display. Protected harmonisation adapters, user inputs and reference
resources were not edited.

## Final regression run

Run from the repository with the local `postgwas` environment. Its Python imports
this repository's `src/postgwas`, verified outside pytest. The following test
selection completed with **675 passed, 1 skipped in 318.17 seconds**:

```bash
R_LIBS_USER=/dev/null python -B -m pytest -q -ra --tb=short -p no:cacheprovider \
  tests/test_validation_consolidation.py \
  tests/test_file_validation_reporting.py \
  tests/test_input_validation.py \
  tests/test_pipeline_validation_reporting.py \
  tests/test_pipeline_validation_presentation.py \
  tests/test_validation_once_migration.py \
  tests/test_filtering.py tests/test_filtering_summary.py \
  tests/test_qc_summary_af.py tests/test_filtering_qc_provenance.py \
  tests/test_flames_handoff.py tests/test_terminal_style.py \
  tests/test_screen_logging.py tests/test_core_helpers.py \
  tests/test_configuration.py tests/test_single_cell_file_validation.py \
  tests/test_ld_annotation.py tests/test_manhattan.py
```

The actual run set `--basetemp` to an invocation-owned temporary directory.
The temporary R setting avoids a broken user-level `data.table` installation
shadowing the environment's working package. No R installation, user library
or persistent environment setting was changed.

| Validation | Expected and observed result |
|---|---|
| Native direct and pipeline LD annotation | Four fixture variants receive the expected two LD-block labels. CHROM, POS, REF, ALT, ES, SE, LP and AF are identical before/after annotation. Malformed BED input fails without a delivered result. |
| Screen and audit inspection | Input filename, build and indexed count occur once in the input-validation card. New output files have independent validation. Hidden mode has empty physical stdout/stderr but retains its transcript and audit. |
| Reuse and failure boundaries | Exact checks run once; stronger contracts and changed direct files revalidate. Pipeline input changes fail. Missing/empty inputs retain their original error policy, including missing cached child dependencies. |
| Reporting edge cases | Same basenames remain distinguishable; conflicting counts and new list-tail values survive; companion failures and screen-limited failures remain visible. |
| Enrichment worker | The same deduplicated gene list reaches the service once; the dispatcher no longer prints a second count. Providers are stubbed in this focused entry-function test. |
| Native Manhattan | Plotting and point-data assertions pass with the environment's working R libraries. Native coding highlighting is the skipped test described below. |
| Shared command boundaries | Registered direct-command progress/audit boundaries, representative pipelines, failure paths, colours, native capture and hidden-screen logging pass. These are not full native runs of every engine. |

No new statistical method was introduced. Scientific validation here means
preservation of the tested data and existing validation boundaries, not universal
validation of every supported analysis. Supporting evidence is in the named
tests and current implementation; no external method change required a new
literature-based decision.

## Outstanding verification and pre-existing issues

1. `tests/test_manhattan.py::test_native_coding_highlight_keeps_unannotated_sites`
   is skipped: a tested split-vep-enabled executable is not configured. The local
   `bcftools plugin -l` lists only `liftover`. Remaining work: provide that
   executable through `POSTGWAS_TEST_SPLIT_VEP_BCFTOOLS` and rerun the native test.
2. `python -B -m pytest -q -ra --tb=short -p no:cacheprovider tests/test_wiki_docs.py`
   reports **18 passed, 1 failed**. The existing `## Locus boundaries` heading in
   `docs/wiki/modules/fine-mapping.md` violates the module-page template's heading
   order. This unrelated page and its valid structural test were not changed.
   Remaining work: reconcile that page's heading hierarchy with the template.
3. The ordinary local R library search still encounters the broken user-level
   `data.table` installation. The temporary test setting is not a permanent fix.
   A separate environment repair or explicit runtime-isolation decision is needed
   before claiming the default local plotting environment is healthy.
4. Full production datasets, all external engines/provider services, independent
   worker runtimes, Linux execution and performance benchmarks were not run.
   The shared checks are reused where registered; this change does not prove
   every opaque engine reads every resource only once. New output validation,
   stricter contracts and analysis reads remain necessary work.

## Compliance handoff

- Checklist: **blocked for unrestricted all-tools/release sign-off** by the
  outstanding checks above; the implemented scope's regression run passes.
- Cumulative task diff reviewed: **yes**, against a pre-edit snapshot so unrelated
  working-tree changes are not confused with this task.
- Hardcoding and duplication audit: **passed for these changes**; presentation
  aliases remain in validated YAML and shared infrastructure owns the display.
- Scientific validation: **completed for the representative invariants above**;
  no new scientific behaviour, no claim of complete all-engine validation.
- Tests run/not run: exact selection, outcomes, skip and remaining risks are above.
- `git diff --check` on the task's modified paths: **passed**.
- Hygiene: only task-owned temporary snapshots and test fixtures are disposable;
  user inputs, references and unrelated changes are preserved.

Accuracy confidence: **95%** for the scoped implementation and observed tests;
the listed native-environment and full-run gaps prevent an all-system guarantee.
