# Logging, progress, and configuration

This module records what it did, which inputs and settings it used, how far it
has progressed, and what remains. User-selectable fine-mapping defaults are
defined once in the schema-validated
`config/defaults/modules/fine_mapping.yaml`; CLI values are suppressed in
`argparse` and only explicit CLI values override matching YAML fields. A run
copies the effective settings to its output directory so results can be audited
without inspecting source code.

Runnable direct and dependency-complete pipeline commands are provided in the
[fine-mapping examples](../../../examples/fine_mapping/README.md).

Both engines use 0.95 as the canonical credible-set coverage and accept an
explicit schema-validated YAML override. FINEMAP also exposes
`--prob-cred-set`; SuSiE reads the resolved common value from the single
configuration document passed to R. Primary and recovery fits always use the
same resolved target.

## Terminal output versus the audit log

The terminal presents eight ordered Rich stages with short outcome counts,
followed by one decision-focused final summary. It reports the engine, locus
success/failure and warning counts, primary and authoritative credible-set
counts, overlap disposition, and the principal result, QC, FLAMES, and log
locations. The output root is printed once and files below it use concise
relative paths. It does not print raw result dictionaries, worker-copy events,
per-locus audit events, resolved paths for every intermediate file, or
`[STAGE]`/`[PROGRESS]` records.

Those details remain available without loss in
`run_metadata/pipeline_summary.log`; machine-readable progress remains in
`run_metadata/pipeline_progress.tsv`, and locus-specific evidence remains in
the QC files and locus logs described below. The global schema-validated
`logging.show_progress` setting controls the live eight-stage display, and
`logging.terminal_label_width` controls its alignment. The CLI introduces no
independent display defaults.

## Progress records

Both engines write `run_metadata/pipeline_progress.tsv` with these fields:

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

## Output layout and retention

The canonical YAML places durable outputs into four meaningful categories and
all working data into one intermediate category:

```text
fine_mapping/
├── results/
│   ├── combined_results/
│   ├── primary_credible_sets/
│   └── diagnostic_plots/
├── quality_control/
│   ├── locus_logs/
│   └── ld_diagnostics/
├── downstream_inputs/
│   └── flames/
├── run_metadata/
└── intermediate_files/
    ├── fitted_models/
    ├── susie_workers/
    ├── locus_chunks/
    ├── prepared_loci/
    ├── summary_statistics/
    ├── finemap_inputs/
    ├── finemap_locus_work/
    ├── selected_finemap_models/
    ├── primary_flames_handoff/
    └── overlap_joint_rerun/
```

`output_layout` in the canonical fine-mapping YAML defines directory placement
and the main engine files; the remaining handoff filenames are also defined in
the schema-validated YAML. The CLI does not carry another set of names. Primary
FLAMES exports, fitted R objects, locus chunks, prepared summary-statistics
tables, chromosome caches, FINEMAP locus work and joint overlap reruns remain
below `intermediate_files/`.

SuSiE workers first publish into the result, QC and metadata paths. After those
outputs validate, packaged `cleanup_successful_workers: true` removes the
duplicated worker tree. Worker or primary-merge failures retain their
intermediate evidence, and setting the policy to `false` explicitly retains
successful workers. Empty
annotation directories are not created by fine-mapping; FLAMES materializes
annotations when it runs.

## FINEMAP audit files

| File | Contents |
|---|---|
| `run_metadata/pipeline_summary.log` | Ordered orchestration stages, counts, paths, warnings, and exceptions. |
| `run_metadata/pipeline_progress.tsv` | Pipeline and locus percentage/remaining records. |
| `run_metadata/run_configuration.json` | Inputs, analysis controls, resource calculation, software identity, and resolved output layout. |
| `quality_control/finemap_input_harmonization_qc.tsv` | Input, invalid, coordinate-mismatch, palindromic, allele-mismatch, sign-flip, and retained counts. |
| `quality_control/finemap_task_generation_skips.tsv` | Loci not submitted and their explicit reason. |
| `quality_control/finemap_locus_status.tsv` | Success/failure status, stable failure code, and detailed reason for every attempted locus, including model-output validation and external-process timeout failures. Timeout rows include stage, configured/elapsed seconds, forced-termination status, and removed partial outputs. |
| `run_metadata/finemap_software_versions.tsv` | PLINK, bgenix, LDstore, FINEMAP, genome build, and resolved coverage. |
| `intermediate_files/finemap_locus_work/<locus>/<locus>_debug.log` | Every per-locus stage, file cleanup, counts, and exception trace. |
| `*.plink.log`, `*.bgenix.log`, `*.ldstore.log`, `*.finemap.log` | Exact external command, stdout, stderr, and exit status. |
| `intermediate_files/finemap_locus_work/<locus>/<locus>_QC.tsv` | Sample/variant counts, configured and effective maximum causal-SNP counts, resolved coverage target, and LD validation statistics. |
| `quality_control/finemap_model_selection_summary.csv` | Every evaluated causal-count model, posterior probability, selection flag, validation status, and model-output failure detail. |
| `intermediate_files/primary_flames_handoff/` | Primary credible sets and manifest retained for overlap provenance. |
| `downstream_inputs/flames/indexfile.txt` | Authoritative post-overlap FLAMES-compatible file/locus/annotation mapping. |

## SuSiE audit files

| File | Contents |
|---|---|
| `run_metadata/pipeline_summary.log` | Python orchestration, splitting, worker collection, merging, validation, cleanup, and exceptions. |
| `run_metadata/pipeline_progress.tsv` | Pipeline and worker percentage/remaining records. |
| `run_metadata/run_configuration.json` | Inputs, analysis controls, resources, software identity, resolved SuSiE configuration, output layout, configuration path, and SHA-256 checksum. |
| `run_metadata/resolved_susie_configuration.json` | Exact schema-validated fine-mapping configuration consumed by every R worker. |
| YAML `summary_statistics_preparation` | Exact-file policy, streaming chunk size, separator, filename templates, and manifest/summary filenames. Directory placement comes only from YAML `output_layout`. |
| `intermediate_files/summary_statistics/chromosome_summary_statistics_manifest.tsv` | Analysed chromosome partitions, row counts, and source summary-statistics provenance. |
| `intermediate_files/summary_statistics/locus_summary_statistics_manifest.tsv` | Exact effective boundary, locus ID, prepared file, variant count, chromosome partition, policy, and source for every locus. |
| `quality_control/summary_statistics_preparation_summary.tsv` | Source rows scanned, analysed chromosomes, prepared loci, variant memberships, cache reuse, and both manifest paths. |
| `quality_control/<sample>_susie_locus_progress.tsv` | Merged R primary/recovery progress for every current-run worker. |
| `run_metadata/<sample>_run_configuration_r.tsv` | Effective R inputs, every resolved YAML control, configuration path, and configuration MD5 checksum. |
| `run_metadata/<sample>_software_versions.tsv` | PLINK, build, and resolved coverage recorded by Python. |
| `run_metadata/<sample>_software_versions_r.tsv` | R and susieR versions, build, and resolved coverage. |
| `run_metadata/<sample>_run_susie.log` | Concatenated current-worker R output. |
| `quality_control/locus_logs/*` | Detailed per-locus R stages, validation, recovery method, and failure reason. |
| `quality_control/<sample>_susie_locus_qc.tsv` | Convergence, recovery, LD mismatch, dense-LD resource estimates, fit-payload serialization/reuse counts, resolved coverage target, and notes by locus. |
| `quality_control/<sample>_susie_failed_loci.tsv` | Failed or skipped loci only. |
| YAML `recovery_audit_filename` | One row for every primary fit, LD validation failure, LD repair/revalidation, and recovery fit, including exact reason, repairability, LD state, model controls, lambda, minimum eigenvalue, and outcome. |
| YAML `sample_size.policy` | Currently `warn`: continue with the configured locus summary while recording material SNP-level NEF heterogeneity. |
| YAML `sample_size.summary_statistic` | Currently `median`, calculated independently within every locus. |
| YAML `sample_size.relative_range_warning_threshold` | Warn when `(maximum NEF - minimum NEF) / median NEF` exceeds this value. The packaged value is `0.05`. |
| `intermediate_files/primary_flames_handoff/` | Primary component-level coverage and mappings retained for overlap provenance. |
| `downstream_inputs/flames/` | Authoritative post-overlap credible sets and FLAMES index. |

Successful SuSiE combined result and credible-set tables, and the FINEMAP FLAMES manifest, include the NEF diagnostics and warning reason. SuSiE combines multiple warnings in the same field—for example, NEF heterogeneity plus a material non-fatal LD/z mismatch—using semicolons. Failed loci cannot have credible-set rows; their exact reasons are retained in `quality_control/*_susie_locus_qc.tsv`, `quality_control/*_susie_failed_loci.tsv`, or `quality_control/finemap_locus_status.tsv` and summarized on screen.

FINEMAP validates every causal-count model's `Post-Pr` header twice: in the
locus worker immediately after FINEMAP exits and again before cross-locus model
selection. A valid value may be exactly zero or one. A missing, malformed,
duplicate, non-finite, or out-of-range value fails the complete affected locus
with `invalid_model_posterior_probability`; it is never converted to zero,
clipped, or omitted from the comparison. Independent valid loci continue and
remain eligible for FLAMES export. The locus status TSV, per-locus log, and
all-model summary retain the stable reason and exact file detail; the final
screen groups affected loci by the stable reason.

FINEMAP external execution is bounded by three schema-validated settings under
`engines.finemap`: `ldstore_timeout_seconds` (21,600 seconds, or 6 hours),
`finemap_timeout_seconds` (43,200 seconds, or 12 hours), and
`termination_grace_seconds` (30 seconds). The LDstore value covers both of its
sequential subprocesses as one per-locus stage budget. These values are
operational defaults, not scientific thresholds, and may be explicitly
overridden in YAML or on the CLI for appropriately benchmarked workloads.

LDstore and FINEMAP run in isolated process groups. When a limit expires, the
worker terminates the complete group, waits for the configured grace period,
and force-kills it if necessary. The per-locus log retains partial stdout and
stderr, elapsed time, and the termination action. The worker removes partial
BCOR/LD or FINEMAP configuration/credible-set outputs, writes a failed locus QC
row with `ldstore_timeout` or `finemap_timeout`, and excludes that locus from
model selection and FLAMES. Other locus workers continue, and the final screen
summary groups timeout counts by their stable reason.

For each scientifically validated SuSiE matrix state, the worker writes the
exact `z`, dense LD matrix, and scalar locus sample size to one read-only
temporary RDS. Fresh child processes still provide the configured hard timeout,
but retries read that same payload and write only a small attempt-control RDS.
The locus log, QC table, and recovery audit report the LD state, payload size,
serialization count, child-fit count, and reuse count. An
original validated matrix and a repaired-and-revalidated matrix are never
represented by the same payload. Temporary payloads are removed after success
or after recovery is exhausted; an unavailable payload fails explicitly rather
than being silently recreated.

A converged fit is not treated as evidence that a credible set exists. If no
sets survive SuSiE's configured purity filter, the module returns
`completed_no_credible_sets`, reports zero credible sets and the number of
converged loci without sets, and returns no `flames_input` path. If only some
converged loci yield sets, `genomic_loci.tsv` and the FLAMES index include only
those loci.

Worker directories are cleared only for pipeline-owned generated artifacts
before a rerun. Merge operations log every copied artifact and validate table
headers, preventing stale or schema-incompatible results from being silently
combined. After primary aggregate validation succeeds, the configured cleanup
policy removes successful worker copies; worker or primary-merge failures
retain their intermediate evidence.

## Dense-LD resource guard

Both engines apply `ld_resource_guard.maximum_variants_per_locus` immediately
before constructing a locus LD matrix. The packaged limit is 30,000 variants,
and the boundary is inclusive: 30,000 may proceed, while 30,001 fails that
locus with `maximum_variants_per_locus_exceeded`. Other primary loci continue.
A failed post-primary joint rerun follows the configured overlap failure policy
and is excluded and reported; its superseded overlapping primary sets are not
silently restored.

The limit never filters or selects SNPs. It stops the entire locus before LD so
that a partial scientific model cannot be mistaken for the requested model.
Every processed locus records its input variant count, configured limit, one
float64 dense-matrix size, engine-specific peak-matrix multiplier, estimated
peak dense-LD memory, reserved worker memory, and resource warning/failure
reason. The matrix size is `variants² × 8 / 1024³` GiB. At 30,000 variants one
matrix is approximately 6.71 GiB; the packaged FINEMAP multiplier of 4 estimates
26.82 GiB peak and the SuSiE multiplier of 5 estimates 33.53 GiB peak. These are
planning estimates, not guaranteed measurements. Worker concurrency is reduced
using the largest estimate, and a warning is retained when the estimate exceeds
reserved worker memory.

Override the canonical YAML value for a planned high-memory analysis or use
`--maximum-variants-per-locus COUNT`. The CLI parser assigns no independent
default. Increasing the limit does not make a large model scientifically more
reliable and can substantially increase memory, text-matrix I/O, and cubic LD
validation cost.

Before workers start, Python makes one streaming pass over the original
summary-statistics table and retains only chromosomes required by the effective
locus boundaries. The configured `validation_scope: analysed_chromosomes`
therefore validates every retained chromosome row and explicitly records the
count outside that scope; rows outside chromosomes used by this analysis do not
affect a locus model. It then loads each retained chromosome once, sorts by `BP`,
and uses inclusive indexed range selection to write one exact table per locus.
The R workers read the small locus manifest and only the exact tables assigned
to their locus chunk; they never load or scan the genome-wide source. If two
primary boundaries overlap, their exact files intentionally contain the shared
variants independently. The post-primary joint rerun reuses the same chromosome
cache and writes a new exact union file with no additional flank.

## Post-primary overlap resolution

`overlap_resolution.policy: rerun_connected_groups` deliberately waits until
the complete primary locus-wise pass has finished. It then forms connected
components from the exact closed boundaries used by successful primary fits.
A converged SuSiE locus remains eligible even when purity filtering retained no
credible set. Non-overlapping primary credible sets are retained unchanged;
each overlapping component is rerun once using the union of its already
expanded boundaries. The configured flank is set to zero for this second pass,
so the union is not expanded again.

The overlap YAML block controls the optional maximum joint span, failure policy,
report filenames, file prefixes, overlap group prefix, and provenance
delimiter. Directory placement is controlled only by `output_layout`.
`maximum_joint_region_kb: null` means
that no universal biological span is assumed; projects may set a positive cap
for their LD reference and compute budget. Under
`joint_failure_policy: exclude_and_report`, a failed or no-credible-set joint
group is excluded from the authoritative final handoff instead of restoring
duplicated primary results.

`overlapping_loci_resolution.tsv` records every primary locus, connected group,
joint boundary, disposition, and reason. `final_combined_credible_sets.tsv`
contains the retained non-overlapping primary rows plus successful joint-rerun
rows. Every row records its analysis round, final selection reason, original
primary locus IDs, original credible-set file paths, primary manifest/status
paths, primary warning/failure reasons, source result, and final copied file.
The returned `flames_input` points to the configured authoritative final
directory, while all primary and joint-run directories remain available for
audit. The configured overlap summary TSV records the primary, overlap, joint
success, no-set, failure/exclusion, and final credible-set counts even when no
authoritative final set can be produced.

## Defaults and overrides

- The canonical FINEMAP baseline is `algorithm: sss`, `n_causal_snps: 5`,
  `n_iter: 100000`, `n_conv_sss: 100`, `prob_conv_sss_tol: 0.001`,
  `n_configs_top: 50000`, `corr_config: 0.95`, `pvalue_snps: 1.0`,
  `cond_pvalue: 5e-8`, `prior_std: 0.05`, and 0.95 credible-set coverage.
  Boolean controls `prior_k`, `force_n_samples`, and `std_effects` default to
  false. These match FINEMAP 1.4.2's documented baseline where applicable.
  `n_causal_snps` is an upper bound, not an exact causal count; for a locus
  containing fewer variants, PostGWAS uses the variant count as the effective
  bound and records both values in the locus log and QC table.
- FINEMAP execution defaults are a 21,600-second total LDstore stage timeout, a
  43,200-second FINEMAP fitting timeout, and a 30-second termination grace
  period. CLI overrides are `--ldstore-timeout-seconds`,
  `--finemap-timeout-seconds`, and
  `--finemap-termination-grace-seconds`; the parser assigns no independent
  values.
- Locus construction, the LP threshold, MHC handling and coordinates, genome
  build, and per-worker memory are common YAML settings consumed by both
  engines. The common `ld_resource_guard`, `validation`, and `runtime` blocks
  resolve dense-LD limits, numerical tolerances, version-probe timeouts, and
  fallback memory. SuSiE's nested `fitting`, `recovery`, `ld_validation`,
  `execution`, `memory`, and `plotting` sections contain every user-selectable
  R control.
- `get_finemap_defaults()` derives its compatibility dictionary from an
  already validated configuration model; it contains no independent numeric
  defaults.
- `susie/defaults.r` contains only protocol invariants and a strict JSON
  loader. Missing resolved controls fail before locus processing; R does not
  supply alternative scientific, recovery, timeout, memory, or plotting
  defaults.
- Help renders the configured-default label in bold green and its YAML value in
  cyan. The displayed value is informational and is not assigned by the parser.
- The selected worker count is derived from requested threads and named RAM
  requirements, then recorded; it is not an unexplained literal in workflow
  code.
- File-format rules, mathematical identities, and required schema names remain
  implementation invariants. They are validated and documented in logs rather
  than treated as tunable scientific thresholds.

FINEMAP option definitions and upstream defaults were checked against the
[FINEMAP 1.4.2 command-line documentation](https://www.christianbenner.com/).
`prior_k` remains disabled unless explicitly requested, and PostGWAS rejects
that request before analysis until its master-file generator supports the
required K-file column. The undocumented historical `--collinear-tol` option
has been removed from the PostGWAS interface.
