# Shared runtime terminal style — validation record

Date: 2026-09-08. Scope: PostGWAS runtime presentation in direct and pipeline
modes. This is not a new full-data GWAS validation or installation audit.

## Behaviour

One shared, schema-validated `logging.terminal_style` now controls runtime
colours and indentation. Active operations are cyan, success green, warnings
yellow, errors red and decisions magenta. Ordinary values and paths use the
terminal foreground. Shared fields align wrapped values; narrow live bars
reserve room for measured counts, percentage and elapsed time. Complete
operation descriptions remain in the durable progress milestones.

Native tool text remains literal and retains its wording and native layout.
Terminal escape sequences cannot override the theme and are removed from
saved transcripts. Unknown native messages remain neutral: the renderer does
not infer severity or validated completion from words in tool output.
`--hide-screen` retains logging; `NO_COLOR` and the YAML `never` setting disable
colour. Cross-stream ordering between concurrent stdout/stderr writers is not
guaranteed; order within each captured stream is preserved.

The change does not alter association calculations, transformations, reference
files, analysis arguments, stage order or output-validation conditions. The
protected harmonisation adapters were compared against the task-start copies
and are unchanged. No existing analysis outputs were deleted or replaced.

## Configuration and implementation

Added under the canonical application YAML:

- `logging.terminal_style.color`: `auto` or `never`, default `auto`;
- `logging.terminal_style.heading_indent`: default `2`;
- `logging.terminal_style.field_indent`: default `6`;
- `logging.terminal_style.indent_step`: default `4`;
- `logging.terminal_style.styles`: the complete semantic-role palette.

The existing `logging.terminal_label_width` is reused. Run overrides are merged
with packaged defaults before schema validation. No CLI defaults or
module-specific theme keys were added. The old module-style override interface
and QC-specific palette were removed, and their callers/tests migrated.

Files changed by this task, relative to the task-start working tree:

- Shared configuration: `src/postgwas/config/defaults/application.yaml`,
  `src/postgwas/config/models/logging.py`.
- Shared presentation/capture: `src/postgwas/core/ui/{__init__,screen,progress}.py`,
  new `src/postgwas/core/ui/terminal_text.py`,
  `src/postgwas/core/screen_logging.py`.
- Shared entry points and messages: `src/postgwas/__main__.py`,
  `src/postgwas/core/checkpointing.py`,
  `src/postgwas/pipeline/{cli,executor}.py`.
- Module CLI message adapters:
  `src/postgwas/modules/{caldera,flames,formatting,gcta_gene,imputation,kpops,ldsc,magma,magmacovar,mixer,pops,single_cell}/cli.py`.
- Remaining module messages: `src/postgwas/modules/enrichment/main.py`,
  `src/postgwas/modules/harmonisation/shared/runtime.py`,
  `src/postgwas/modules/manhattan/service.py`,
  `src/postgwas/modules/qc_summary/{reporting,service}.py`.
- Tests: new `tests/test_terminal_style.py`, `tests/test_core_helpers.py`,
  `tests/test_qc_summary_af.py`, `tests/test_cli_help.py`.
- Rules/documentation: `AGENTS.md`,
  `docs/wiki/core/logging-and-reproducibility.md`, and this record.

Scientific decisions: none introduced. Terminal implementation was checked
against installed Rich APIs and the official [console documentation](https://rich.readthedocs.io/en/stable/console.html),
[live-rendering documentation](https://rich.readthedocs.io/en/stable/live.html)
and [style documentation](https://rich.readthedocs.io/en/stable/style.html).

## Validation

Environment: macOS, Python 3.10 in the local `postgwas` environment. Commands
below used `/Users/JJOHN41/miniconda3/envs/postgwas/bin/python` as `python`.
Results overlap and must not be added together as a unique test count.

Broad focused regression: **449 passed, 1 skipped**.

```bash
python -B -m pytest -q -p no:cacheprovider -ra \
  tests/test_terminal_style.py tests/test_screen_logging.py tests/test_core_helpers.py \
  tests/test_qc_summary_af.py tests/test_harmonisation_pvalue_inference_warning.py \
  tests/test_manhattan.py tests/test_configuration.py tests/test_checkpointing.py \
  tests/test_direct_execution.py tests/test_pipeline_validation_presentation.py \
  tests/test_pipeline_validation_reporting.py tests/test_pipeline_pilot_regressions.py \
  tests/test_installation_script.py \
  tests/test_cli_help.py::test_pipeline_error_renders_bracketed_external_tool_paths_literally \
  tests/test_cli_help.py::test_pipeline_gcta_missing_external_inputs_stop_before_execution \
  tests/test_cli_help.py::test_finemap_pipeline_bare_selection_stops_before_analysis
```

Latest shared-presentation regression after the final narrow-bar refinement:
**175 passed**.

```bash
python -B -m pytest -q -p no:cacheprovider \
  tests/test_terminal_style.py tests/test_core_helpers.py tests/test_screen_logging.py
```

Dedicated colour checks with inherited `NO_COLOR` removed: **94 passed**.

```bash
env -u NO_COLOR python -B -m pytest -q -p no:cacheprovider tests/test_terminal_style.py
```

Representative fixtures and expected/observed results:

- All 21 registered scientific direct-command boundaries, successful and failed
  analysis stubs: shared presentation and screen recording observed; failures
  did not become completed runs. These are dispatcher tests, not 21 full GWAS runs.
- Registered pipeline stages with stub runners/checkpoints: execution order,
  success and second-stage failure matched expected progress accounting.
- Real CLI error adapters with injected configuration errors: actionable,
  literal paths including square brackets; shared error styling and stderr
  behaviour matched expectations.
- Narrow/wide terminal fixtures (60, 80, 120 and 180 columns): shared field
  continuation alignment and visible progress counts, percentage and elapsed
  time matched expectations. Terminal and forced-colour output were inspected.
- Subprocess and nested-progress fixtures: visible, redirected and hidden
  output; UTF-8 split across chunks; ANSI/OSC controls; CRLF; malformed control
  sequences; failure and drain handling. Printable diagnostics were retained,
  saved logs had no terminal controls, and hidden output still logged.
- Configuration fixtures: valid overrides, context restoration, unknown roles,
  incomplete palettes and invalid styles produced the expected results.
- Syntax parsing passed for all 31 task-changed Python files. `git diff --check`
  passed. Removed palette/override names have no remaining implementation
  callers. No repository-configured lint/type-check command was identified.

## Outstanding checks and operational limitations

An earlier broader CLI-help audit returned **219 passed, 2 failed**. One failure
was the old expected pipeline error prefix; it was updated for the requested
shared error symbol without weakening the literal-path assertion and passed in
the final focused run. The other remains:

- `tests/test_cli_help.py::test_readme_prose_avoids_unnecessary_scientific_wording`
  rejects the pre-existing wording in `README.md` around line 240:
  “requirements fail before chromosome work or scientific outputs begin.”
  This task did not edit README. Correcting that unrelated prose and rerunning
  the test remain outstanding; it is not a runtime-style failure.

The skipped test is
`tests/test_manhattan.py::test_native_coding_highlight_keeps_unannotated_sites`.
Installed bcftools reports version 1.23.1, but `bcftools +split-vep -h` fails to
load the plugin in this environment. Testing that native coding-highlight path
requires a verified split-vep-enabled executable and
`POSTGWAS_TEST_SPLIT_VEP_BCFTOOLS`. No plugin installation was attempted for this
presentation-only task.

Full-scale runs of every external analysis tool, Linux execution, all terminal
emulators, and end-to-end performance benchmarking were not performed. Shared
boundary tests establish presentation behaviour, not scientific validity of
every external analysis. Existing checkpoint policy fingerprints package and
configuration content: installing these changes can therefore trigger the
existing warning-and-restart policy on the next run. That policy was not changed.

The working tree already contained substantial unrelated changes. Review used
task-start copies to isolate this task's cumulative changes. Existing result
directories and screen logs were preserved, including diagnostic lines appended
by CLI validation tests. Only task-created temporary audit copies and the
temporary display demonstration are cleanup targets.

## Compliance handoff

- Checklist: **blocked for an all-checks-passed claim**, by the outstanding
  README check and unavailable optional native plugin; focused presentation
  regressions passed.
- Cumulative diff reviewed: **yes**, against the task-start working tree.
- Hardcoding and duplication audit: **passed** for the task changes; runtime
  defaults are in YAML and presentation is shared.
- Scientific validation: **not applicable to new scientific behaviour**;
  no scientific algorithm was changed and the relevant existing regression
  suites passed subject to the explicit native-test skip.
- Tests run: commands/results above.
- Tests not run and outstanding issues: explicitly listed above. No claim of
  a clean full-repository test run or all-platform production validation.
