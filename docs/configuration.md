# Configuration architecture

PostGWAS resolves configuration in this order:

1. Packaged per-module defaults.
2. An optional packaged profile.
3. User YAML, including optional relative `include` files.
4. Explicit CLI overrides.

Module help uses the shared Rich formatter for defaults: the default label is
shown in bold green and the resolved value in cyan. This is presentation only;
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

Resume is enabled by the single canonical `run.resume` default.
Resume-aware direct modules and pipeline steps therefore reuse only outputs that
pass their completion checks; `--no-resume` disables continuation for one run.
After a formatter step completes, PostGWAS atomically records the input VCF
checksum, resolved formatter-configuration checksum, output checksums, and
reusable result contract. A repeated pipeline validates that completion manifest
and continues with the remaining steps. Changed, missing, empty, or
checksum-mismatched inputs and outputs are rejected rather than silently reused.
When a completed formatter run directory is copied, resume resolves artifacts
from the current output directory, requires their relative paths to match the
current configured filenames, validates their checksums there, and rewrites the
copied manifest with current absolute paths. It never returns an artifact from
the manifest's former output directory.
Use `--no-resume` or set `run.resume: false` to force collision-protecting
non-resume behaviour. Set `run.overwrite: true` or pass `--overwrite` to rerun a
step and replace its outputs; overwrite takes precedence over resume.
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
validated scientific policy registry. There is no second defaults file inside
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
# Complete explanations for every harmonisation policy
postgwas config export --module harmonisation --style full > harmonisation.yaml

# Short comments
postgwas config export --module harmonisation --style minimal > harmonisation.yaml

# Key-value pairs only; the global shorthand is equivalent
postgwas --config --module harmonisation --style values > harmonisation.yaml

# Fine-mapping plus every required preceding module, in execution order
postgwas config export --pipeline finemap > finemap_pipeline.yaml

# Combine several final analysis targets without duplicating prerequisites
postgwas config export --pipeline finemap magma > analysis_pipeline.yaml
```

Use `--output PATH` instead of shell redirection when preferred. Pipeline
export can also resolve values from an existing run file with `--run-config PATH`.
The dependency planner determines all preceding modules; no profile selection
is required. Method-dispatched modules also retain non-executable configuration
they consume. For example, an LDSC cell-type pipeline exports `formatting`, the
supporting `ldsc` munging settings, and `single_cell`, while only `formatting`
and `single_cell` appear in the executable `pipeline.modules` list.
