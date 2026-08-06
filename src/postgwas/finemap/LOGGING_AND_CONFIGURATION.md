# Logging, progress, and configuration

This module records what it did, which inputs and settings it used, how far it
has progressed, and what remains. Scientific, numerical, resource, timeout,
and plotting controls are named in `defaults.py` (Python orchestration and
FINEMAP) or `susie/defaults.r` (SuSiE/R). A run copies the effective settings
to its output directory so results can be audited without inspecting source
code.

The credible-set target is intentionally different from other defaults: it is
fixed at 0.95 in both defaults files and runtime validation rejects a different
target. This prevents a caller from silently producing a non-95% set.

## Progress records

Both engines write `pipeline_progress.tsv` with these fields:

| Field | Meaning |
|---|---|
| `timestamp_utc` | UTC time at which the event was recorded. |
| `scope` | Unit being tracked, such as the whole pipeline, loci, or workers. |
| `stage` | Named operation being started, completed, skipped, or failed. |
| `status` | Outcome at the time of the event. |
| `completed` | Completed units within this scope. |
| `total` | Total units within this scope. |
| `percentage` | `100 * completed / total`, rounded only for display. |
| `remaining` | Exact integer `total - completed`. |
| `message` | Locus, chunk, output, counts, or failure reason. |

The same event is written to the human-readable pipeline log as a
`[PROGRESS]` entry. A failed run retains the last completed percentage and
remaining count and adds a final `failed` event; it does not claim 100%.

FINEMAP pipeline stages are initialization, dependency validation, input
preparation, task generation, locus execution, output formatting, and
completion. SuSiE uses initialization, dependency validation, locus
preparation, locus splitting, worker execution, merge/validation, and
completion. Separate scope rows report per-locus or per-worker progress.

## FINEMAP audit files

| File | Contents |
|---|---|
| `pipeline_summary.log` | Ordered orchestration stages, counts, paths, warnings, and exceptions. |
| `pipeline_progress.tsv` | Pipeline and locus percentage/remaining records. |
| `run_configuration.json` | Inputs, analysis controls, resource calculation, software identity, and every named Python default. |
| `input_harmonization_qc.tsv` | Input, invalid, coordinate-mismatch, palindromic, allele-mismatch, sign-flip, and retained counts. |
| `task_generation_skips.tsv` | Loci not submitted and their explicit reason. |
| `finemap_locus_status.tsv` | Success/failure status and reason for every attempted locus. |
| `software_versions.tsv` | PLINK, bgenix, LDstore, FINEMAP, genome build, and fixed coverage. |
| `loci_data/<locus>/<locus>_debug.log` | Every per-locus stage, file cleanup, counts, and exception trace. |
| `*.plink.log`, `*.bgenix.log`, `*.ldstore.log`, `*.finemap.log` | Exact external command, stdout, stderr, and exit status. |
| `loci_data/<locus>/<locus>_QC.tsv` | Sample/variant counts, fixed target, and LD validation statistics. |
| `finemap_best_results/finemap_all_models_summary.csv` | Every evaluated causal-count model, posterior probability, and selected model. |
| `flames_input/finemap_FLAMES_manifest.tsv` | Source model, component, fixed target, achieved coverage, variant count, build, and allele definition. |
| `flames_input/indexfile.txt` | FLAMES-compatible file/locus/annotation mapping. |

## SuSiE audit files

| File | Contents |
|---|---|
| `pipeline_summary.log` | Python orchestration, splitting, worker collection, merging, validation, and exceptions. |
| `pipeline_progress.tsv` | Pipeline and worker percentage/remaining records. |
| `run_configuration.json` | Inputs, analysis controls, resources, software identity, and every named Python default. |
| `<sample>_SuSiE_locus_progress.tsv` | Merged R primary/recovery progress for every current-run worker. |
| `<sample>_run_configuration_r.tsv` | Effective R arguments plus every named SuSiE default. |
| `<sample>_software_versions.tsv` | PLINK, build, and fixed coverage recorded by Python. |
| `<sample>_software_versions_r.tsv` | R and susieR versions, build, and fixed coverage. |
| `<sample>_run_susie.ALL.log` | Concatenated current-worker R output. |
| `logs/*` | Detailed per-locus R stages, validation, recovery method, and failure reason. |
| `<sample>_SuSiE_QC_summary.tsv` | Convergence, recovery, LD mismatch, fixed target, and notes by locus. |
| `<sample>_SuSiE_failed_loci.tsv` | Failed or skipped loci only. |
| `flames_input/*_manifest.tsv` and `indexfile.txt` | Component-level coverage and FLAMES-compatible mappings. |

Worker directories are cleared only for pipeline-owned generated artifacts
before a rerun. Merge operations log every copied artifact and validate table
headers, preventing stale or schema-incompatible results from being silently
combined.

## Defaults and overrides

- `get_finemap_defaults()` returns a copy of every Python default for logging,
  inspection, and task propagation.
- `get_susie_defaults()` returns all R scientific, retry, resource, timeout,
  memory, MHC, and plotting defaults.
- Functions expose applicable controls as named default arguments. CLI values
  remain supported and are written to the run configuration.
- The selected worker count is derived from requested threads and named RAM
  requirements, then recorded; it is not an unexplained literal in workflow
  code.
- File-format rules, mathematical identities, and required schema names remain
  implementation invariants. They are validated and documented in logs rather
  than treated as tunable scientific thresholds.

