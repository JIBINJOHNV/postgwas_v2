# Configuration

PostGWAS keeps its defaults in schema-validated YAML. For a simple run, use the
public CLI options for your chosen module or pipeline. Use a run configuration
when you want to reuse settings, specify advanced policies, or provide values
that have no public CLI override. Export a template rather than writing one from
memory; the existence of a configuration key does not mean every command exposes
it on the CLI.

## Which configuration file do I need?

| File or interface | Purpose | What it does not do |
|---|---|---|
| Harmonisation sample sheet | Identifies each raw dataset and maps its actual columns | Does not replace the analysis/resource YAML |
| Module run configuration | Supplies supported settings and input/reference paths for one direct analysis | Does not cause the global planner to create missing upstream analyses |
| Pipeline run configuration | Supplies settings for the selected targets and their prerequisites | Does not download references or execute simply because it was exported |
| CLI overrides | Changes the corresponding supported setting for this invocation | Does not make an incompatible file or undeclared method valid |

Most direct commands accept `--run-config`. Direct `manhattan` and
`pathway_enrichment` currently do not; pipeline Manhattan does. Enrichment's
registered `providers` setting is not an implemented provider-selection
interface. Follow the exact command's help rather than assuming that an exported
YAML can drive every direct command. The pipeline requires target selection
through `--modules` or the documented workflow switches.

## Create a configuration file

Export the settings for one module:

```console
postgwas config export \
  --module filtering \
  --style full \
  --output filtering.yaml
```

Limit a formatter module export to one or more downstream formats:

```console
postgwas config export \
  --module formatting \
  --format ldsc \
  --style minimal \
  --output formatting.yaml
```

The resulting module-only YAML contains shared formatter settings and the LDSC
sections, including `exports.ldsc`, but omits other target schemas. Omit
`--format` for the complete formatter configuration. This selector is valid
only with `--module formatting`.

Export a complete pipeline configuration, including required preceding stages:

```console
postgwas config export \
  --pipeline finemap magma \
  --output analysis.yaml
```

Edit the exported file to add your resource paths and study-specific choices.
The [Configuration Defaults](../reference/configuration-defaults.md) page shows
the packaged values for the current checkout.

## Which value is used?

When the same setting is supplied more than once, PostGWAS applies this order:

1. Packaged defaults.
2. Values in your YAML file, including any files it pulls in with `include:`.
3. Explicit command-line overrides.

The corresponding supported command-line value therefore has the highest
priority. Hardware-dependent values such as `threads: auto` and
`memory_gb: auto` are resolved for the machine running the analysis.

## Validate before analysis

Validate a run file without starting the pipeline:

```console
postgwas config validate --config analysis.yaml
```

Unknown keys, invalid values, malformed YAML, missing included files, and
include cycles are rejected. Correct these errors before launching a long run.

## Inspect the effective settings

Inspect the packaged defaults combined with the supplied YAML:

```console
postgwas config show --config analysis.yaml
```

Inspect one module only:

```console
postgwas config show \
  --config analysis.yaml \
  --module fine_mapping
```

Save the effective configuration with the results:

```console
postgwas config show \
  --config analysis.yaml \
  --output results/run_metadata/resolved_config.yaml
```

This inspection command does not replay extra analysis CLI overrides from a
different command. Keep the actual run's resolved configuration and exact
command with the logs, input versions, and reference versions; those record
the settings consumed during execution.

## Screen display and transcript

`logging.show_screen` defaults to `true`. All scientific module and pipeline
commands accept `--hide-screen` as a one-way command-line override. PostGWAS
always appends the complete standard-output and standard-error stream to the
relative `logging.screen_log_file` under `run.output_directory`; its packaged
location is `run_metadata/screen.log`. `--hide-screen` suppresses the terminal
copy only, so the transcript and canonical module logs are still written. Set
`logging.show_screen: true` in a supported run configuration to enable display
when that configuration has disabled it.

`logging.progress_refresh_seconds` controls how often PostGWAS samples
append-only native outputs when an external tool exposes measurable work units.
It defaults to `1.0`; each refresh reads only bytes appended since the preceding
sample. A module may count complete appended result rows or consume an upstream
native work-unit counter, according to when that tool writes its results.

## Resume and overwrite

The same schema-validated policy applies to every scientific direct command and
every pipeline stage:

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

- A completed stage is reused only after its identity, resolved parameters,
  tracked inputs, software, upstream dependencies, declared outputs, SHA-256
  checksums, and scientific completion state validate.
- An unchanged input may reuse its recorded SHA-256 only when its canonical
  path, device, inode, size, modification time, and change time all match the
  checkpoint. Any mismatch forces a fresh content hash, and an input change
  during execution prevents checkpoint publication.
- A real partial run resumes completed pipeline stages or validated internal
  module stages and reruns the incomplete stage from its earliest safe boundary.
- Changed parameters, software, tracked inputs, or upstream checkpoints produce
  a terminal and audit-log warning, then automatically invalidate and rerun the
  affected boundary and its downstream stages.
- Automatic cleanup is limited to checksum-matching PostGWAS-owned regular files
  inside the configured output root. Modified outputs, symlinks, unknown files,
  paths outside that root, and missing tracked inputs are preserved and cause an
  actionable error instead of unsafe replacement.
- Set `run.resume: false` in YAML to disable reuse without requesting output
  replacement. `run.overwrite: true` or `--overwrite` starts the analysis from
  its first step and replaces files created by the previous run. The CLI
  `--resume` and `--overwrite` options are mutually exclusive.

Resume and overwrite are execution controls and do not enter the scientific
checkpoint content digest or get restored from an earlier pipeline stage. A
successful forced overwrite can therefore be reused by the next normal run;
changes to content-determining parameters, inputs, resources, software, or
upstream checkpoints still invalidate reuse.

Do not manually combine outputs from different configurations or studies.

## Related pages

- [Quick Start](../getting-started/quick-start.md)
- [Running Modules Independently](running-modules-independently.md)
- [Pipeline Workflow](pipeline-workflow.md)
- [Input and Output Contracts](input-output-contracts.md)
