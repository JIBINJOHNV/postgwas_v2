# Configuration

PostGWAS uses YAML configuration files for analysis settings, reference
resources, output behavior, and compute limits. Start from an exported template
instead of writing a configuration file from memory.

## Create a configuration file

Export the settings for one module:

```console
postgwas config export \
  --module filtering \
  --style full \
  --output filtering.yaml
```

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
2. An optional packaged profile.
3. Values in your YAML file.
4. Explicit command-line overrides.

The command-line value therefore has the highest priority. Hardware-dependent
values such as `threads: auto` and `memory_gb: auto` are resolved for the
machine running the analysis.

## Validate before analysis

Validate a run file without starting the pipeline:

```console
postgwas config validate --config analysis.yaml
```

Unknown keys, invalid values, malformed YAML, missing included files, and
include cycles are rejected. Correct these errors before launching a long run.

## Inspect the effective settings

Show the values that will be used after defaults and overrides are combined:

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

Keep this resolved file with the command, logs, input versions, and reference
resource versions needed to reproduce the analysis.

## Resume and overwrite

- `run.resume: true` allows supported stages to reuse outputs that pass their
  completion and provenance checks.
- `--no-resume` disables reuse for one command.
- `run.overwrite: true` or `--overwrite` permits a stage to replace its prior
  outputs. Overwrite takes priority over resume.

Do not combine outputs from different configurations or studies. If inputs,
resources, or scientifically important settings change, use a new output
directory or explicitly rerun the affected stages.
