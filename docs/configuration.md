# Configuration architecture

PostGWAS resolves configuration in this order:

1. Packaged per-module defaults.
2. An optional packaged profile.
3. User YAML, including optional relative `include` files.
4. Explicit CLI overrides.

Module help uses the shared Rich formatter for defaults: the default label is
shown in bold green and the resolved value in green. This is presentation only;
configuration-connected CLI parsers keep configurable defaults suppressed so
YAML remains the runtime source of truth. This is the required convention for
every module migrated to the canonical configuration service; legacy parsers
must not be converted by suppressing defaults alone because their services must
first receive the resolved YAML values.

Pipeline help exposes shared workflow inputs and the external resources and
settings used by the selected analyses. Artifacts created by preceding pipeline
steps remain hidden because the pipeline supplies them automatically. Standalone
module help additionally shows the module-specific input artifacts that the user
must provide directly.

Terminal display is enabled by the canonical `logging.show_screen: true`
default. Every scientific module and pipeline accepts `--hide-screen` as a
one-way command-line override for that run. Visible and hidden runs append the
complete standard-output and standard-error stream to
`<run.output_directory>/<logging.screen_log_file>`. The packaged transcript
path is `run_metadata/screen.log`, and schema validation requires it to remain
relative to the output directory. `--hide-screen` changes display only and does
not disable the transcript or the module's canonical scientific logs. To
enable display when YAML has disabled it, set `logging.show_screen: true` in the
run configuration.
Measured external-tool progress is refreshed from append-only native output at
the schema-validated `logging.progress_refresh_seconds` interval (default: one
second). The monitor reads only newly appended bytes: it can count completed
result rows when the upstream tool writes them during computation, or consume
the tool's native measured counter when results are written only at the end.

Resume is enabled globally by `run.resume: true`. Every scientific direct
command and every pipeline stage uses `run.resume_policy`; there is no separate
module-specific default. PostGWAS atomically records the resolved-configuration
digest, tracked input and output SHA-256 checksums, software identity,
scientific completion status, and upstream checkpoint dependencies. A completed
boundary is reused only when all of that evidence still matches.

For large immutable references, the checkpoint also records the file device,
inode, size, modification time, and change time observed with each SHA-256.
PostGWAS reuses the recorded content hash only while that complete filesystem
identity is unchanged. Any identity change triggers a fresh SHA-256 calculation;
an input that changes during an active execution prevents checkpoint creation.
This avoids repeatedly reading multi-gigabyte reference panels without weakening
the tracked content fingerprint.

The packaged policy is:

```yaml
run:
  resume: true
  overwrite: false
  resume_policy:
    checkpoint_validation: sha256
    checkpoint_directory: run_metadata/checkpoints
    direct_manifest: '{command}_direct.yaml'
    pipeline_stage_manifest: '{stage_number}_{module}.yaml'
    audit_log: checkpoint_events.log
    partial_results: resume_validated_stages
    changed_parameters: warn_and_restart
    changed_inputs: warn_and_restart
    unvalidated_outputs: warn_and_restart
```

On a real partial run, completed pipeline stages and validated native-module
sub-stages are resumed. The incomplete stage restarts from its earliest safe
boundary; a direct module without a validated internal sub-stage restarts the
whole module. If resolved parameters, software, tracked inputs, or an upstream
checkpoint changed, PostGWAS prints and logs a warning, invalidates affected
downstream checkpoints, removes only checksum-matching PostGWAS-owned files
inside the run output directory, and reruns automatically. A missing or
externally modified tracked input, a modified output, a symlink, an unknown
artifact, or a path outside the output root is not automatically replaced
because ownership and scientific integrity cannot be proved.

Set `run.resume: false` in YAML to disable checkpoint reuse while retaining
normal collision protection. `run.overwrite: true` or `--overwrite` starts the
analysis from its first step and replaces files created by the previous run.
The command-line `--resume` and `--overwrite` options are mutually exclusive;
an overwrite setting in YAML takes precedence over resume. Automatic policy
restart does not silently change the resolved overwrite setting. Decisions and
warnings are appended to the configured checkpoint audit log.

Resume and overwrite are execution controls, not result-defining parameters, so
they are excluded from the checkpoint content digest and are never restored from
an earlier pipeline stage. This means a successful `--overwrite` run is reused
normally by the next default-resume invocation. Scientific/module parameters,
resources, tracked inputs, software identity, and upstream checkpoints remain
content determining and invalidate reuse when changed.

Formatter and native module manifests can record more detailed scientific
contracts than the global command boundary. For example, a formatter checkpoint
also records content-determining external inputs such as an active LDSC
merge-alleles table and its reusable result contract.
When a completed formatter run directory is copied, resume resolves artifacts
from the current output directory, requires their relative paths to match the
current configured filenames, validates their checksums there, and rewrites the
copied manifest with current absolute paths. It never returns an artifact from
the manifest's former output directory.
For MAGMA and GCTA formatter outputs created before completion manifests were
introduced, PostGWAS requires a completed canonical log, the same recorded VCF
path, an unchanged resolved formatter configuration, older-than-output VCF
modification time, exact configured headers, complete records, and nonempty
tables. Outputs created from one shared validation set, such as MAGMA SNP-location
and p-value tables, must also have equal variant counts. PostGWAS then establishes
the checksummed manifest used by subsequent runs.

A completed MAGMA analysis from before the gene-level multiple-testing report
was introduced can also be resumed without repeating MAGMA. PostGWAS requires
the canonical MAGMA `COMPLETED` record and every previously required published
artifact; the configured corrected-gene report must be the only missing output.
It then validates the existing MAGMA gene table and derives the configured
Bonferroni/FDR columns atomically. Any other missing output remains an incomplete
run and is not repaired automatically.

Values live under `src/postgwas/config/defaults/`. Pydantic models under
`src/postgwas/config/models/` define types and validation but do not contain
scientific defaults. Filtering and harmonisation are connected to the resolver.
During migration, each remaining scientific module must be changed to receive
a validated configuration object and stop defining fallback values internally.

Harmonisation has one packaged source of truth:
`src/postgwas/config/defaults/modules/harmonisation.yaml`. Its `module` section
contains module-level defaults, while the remaining sections contain the
validated scientific policy registry and schema-checked export presentation
metadata. There is no second defaults file inside
the harmonisation implementation package.

Validate a complete file:

```console
postgwas config validate --config examples/configs/pipeline_standard.yaml
```

Inspect one resolved module:

```console
postgwas config show --module filtering \
  --config examples/configs/modules/filtering.yaml
```

Write a reproducible resolved configuration:

```console
postgwas config show --config pipeline.yaml \
  --output results/run_metadata/resolved_config.yaml
```

Unknown keys, invalid ranges, malformed YAML, and include cycles are rejected.
Machine-specific reference paths belong in user configuration, not in source
code or packaged defaults.

Resolved metadata written by a scientific module is deliberately scoped. It
contains the resolved run, execution, and logging settings; only the global
resources consulted by that module; and only the module configuration consumed
by that step. A related module is retained only when its schema is read during
the analysis, such as the formatter contract validated by GCTA. These scoped
documents remain reloadable because omitted sections are supplied by the same
canonical packaged defaults. `postgwas config show --output`, which is an
application-level inspection command rather than run metadata, continues to
write the complete resolved application configuration.

Formatter schemas follow the same rule. VCF query expressions and every
downstream column mapping are stored only in
`src/postgwas/config/defaults/modules/formatting.yaml`. A run configuration may
override them under `modules.formatting`; the typed model rejects unknown
canonical columns, incomplete target schemas, and invalid transformations.
Before VCF extraction, formatter preflight renders every selected output path,
including named and chromosome-partitioned outputs, and rejects any duplicate
destination. Per-step formatter metadata retains only the selected target
schemas—for example, a MAGMA pipeline records `exports.magma` but omits the
GCTA, fine-mapping, imputation, LDSC, and MiXeR export definitions. The metadata
field selection is itself validated from the same canonical formatter YAML
rather than being duplicated in module code.

The formatter's optional custom CLI table follows the same architecture without
requiring users to edit YAML. `--custom-output` and field header options become
explicit overrides of an otherwise inactive `custom_output` section. The
packaged configuration owns each field's canonical source, transformation, and
validity rule; the CLI supplies only the output filename, requested roles, and
header names. Active custom settings are included in scoped resolved metadata
and completion-manifest checks, while inactive runs omit that section.

## Export configuration templates

Every export is generated from the same validated defaults. Harmonisation
policies follow analysis order, and sequences such as chromosome sets are kept
on one line.

```console
# Common harmonisation settings, including run controls and post-merge QC rules
postgwas config export --module harmonisation --scope common \
  --output harmonisation.yaml

# Complete explanations for every harmonisation policy
postgwas config export --module harmonisation --style full > harmonisation.yaml

# Short comments
postgwas config export --module harmonisation --style minimal > harmonisation.yaml

# Key-value pairs only; the global shorthand is equivalent
postgwas --config --module harmonisation --style values > harmonisation.yaml

# LDSC formatter settings without unrelated formatter target schemas
postgwas config export --module formatting --format ldsc \
  --style minimal --output formatting.yaml

# Fine-mapping plus every required preceding module, in execution order
postgwas config export --pipeline finemap > finemap_pipeline.yaml

# Combine several final analysis targets without duplicating prerequisites
postgwas config export --pipeline finemap magma > analysis_pipeline.yaml
```

Use `--output PATH` instead of shell redirection when preferred. Pipeline
export can also resolve values from an existing run file with `--run-config PATH`.

### A shorter harmonisation configuration

Use `--scope common` for a practical starting file. Set `resources.root` and
`run.output_directory`, then run:

```console
postgwas harmonisation --sample-sheet studies.csv --run-config harmonisation.yaml
```

The common export is a **full-run-shaped** file: global `run`, `execution`,
`resources`, and `logging` settings are followed by `modules.harmonisation`
and the supporting `modules.qc_summary.rules`. Dataset paths and source-column
mappings remain in the sample sheet. In contrast, `--scope all` (the existing
default) exports the **module-only** configuration, with `policies` at its root.
Both shapes are accepted by harmonisation; existing key names are unchanged.

Within harmonisation, the order is reference selection → dataset/chromosome
policies → post-merge checks → advanced file, column and adapter definitions.
Policy groups name the processing-log steps that consume them; settings used at
several steps appear only once. Comments distinguish variant rejection during
harmonisation from the report-only post-merge QC assessment. The palindromic
ambiguity interval is not a blanket filter on every variant's EAF.

`--scope` controls **which settings** appear; `--style full|minimal|values`
controls **how much explanation** appears. `minimal` is the default style.
Common scope is currently supported only by `--module harmonisation`; pipeline
and other module exports remain complete.

To shorten an existing configuration without losing its resolved custom values:

```console
postgwas config export --module harmonisation --scope common \
  --run-config existing.yaml --output harmonisation.yaml
```

Every non-default value is retained, even an advanced setting not in the common
selection, or a global/supporting-module override. Includes are resolved first;
the exported file does not depend on the original include files. Consequently,
a heavily customized common export may be longer than the default template.
Omitted settings still inherit the **installed version's** defaults: archive the
run's resolved configuration and software version for reproducibility. Exported
compute budgets are resolved for this computer; use `auto` for `execution.threads`
and `execution.memory_gb` when you want them resolved again on another computer.

To inspect every harmonisation setting, including advanced definitions and full
policy explanations, export `--scope all --style full`. To archive every resolved
run setting (including global and supporting modules), use
`postgwas config show --config harmonisation.yaml --output resolved.yaml`.
This restructuring changes presentation and export selection, not analysis
defaults, allele handling, filtering, or transformations.

### Other module and pipeline exports

For `--module formatting`, `--format FORMAT [FORMAT ...]` writes shared
formatter settings plus only the selected target-specific sections. Omit it to
retain the complete formatter configuration. The option is rejected with
`--pipeline` and with modules other than `formatting`.
The dependency planner determines all preceding modules; no profile selection
is required. Method-dispatched modules also retain non-executable configuration
they consume. For example, an LDSC cell-type pipeline exports `formatting`, the
supporting `ldsc` munging settings, and `single_cell`, while only `formatting`
and `single_cell` appear in the executable `pipeline.modules` list.
