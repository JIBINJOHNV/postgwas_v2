# Logging and reproducibility

A reproducible PostGWAS result must identify the input data, effective
configuration, software commands, reference resources, transformations, and
completion status that produced it. Terminal progress alone is not an audit
record.

## Two output audiences

PostGWAS separates concise terminal progress from detailed file logging:

- terminal output shows ordered stages, important warnings, failures, and the
  final status needed to monitor a run;
- the append-only `logging.screen_log_file` transcript preserves the complete
  standard-output and standard-error stream beneath the run output directory;
- canonical logs retain detailed parameters, decisions, counts, commands,
  outputs, and failure context.

Terminal display is on by default (`logging.show_screen: true`). Every
scientific module and pipeline accepts `--hide-screen` to suppress terminal
display and `--show-screen` to override a hidden YAML setting. Either choice
continues writing the screen transcript; it does not replace or disable the
structured canonical logs.

On an interactive terminal, progress uses one transient live region. A
scientifically detailed stage display temporarily replaces its enclosing
module or pipeline bar, then returns control when the detailed stages finish;
nested bars never refresh independently. When output is redirected, captured,
or hidden, PostGWAS writes one deterministic `Started`, `Completed`, or `Failed`
milestone per event instead of appending terminal-redraw frames to the durable
screen log.

When an external tool exposes measured work through incrementally written rows
or a native compute counter, PostGWAS adds a nested counter while that tool is
active. The live and durable displays report completed units, the tool-declared
total, percentage, and elapsed time. Intermediate refreshes never reach 100%;
completion is shown only after the external command and its required result
validation succeed. The append-only monitor interval comes from
`logging.progress_refresh_seconds`.

Formatter and MAGMA key-value fields use the shared
`logging.terminal_label_width` setting, so labels and colons remain aligned
across MAGMA stage outcomes and both completion summaries. Wrapped values
continue beneath the value column.

## Resolved configuration

Configuration is resolved once using the documented precedence order. A
scientific run should record the run, execution, logging, resource, and module
values it actually consumed. This is more reliable than preserving only the
user-provided YAML because packaged defaults and explicit CLI overrides also
affect execution.

Inspect or save resolved values with the configuration command:

```console
postgwas config show \
  --config pipeline.yaml \
  --output results/run_metadata/resolved_config.yaml
```

## Minimum provenance for a completed stage

The run record should make the following information recoverable:

- dataset identity and input paths;
- input metadata and validation results;
- effective configuration and scientific policies;
- external commands and software versions;
- reference files, genome build, and population context;
- row or variant counts before and after transformations;
- exclusions, inferred or transformed values, and their reasons;
- generated artifacts and their validation status;
- per-stage runtime and configured compute resources;
- warnings, failures, and final completion status.

Do not interpret the presence of an output file as proof of success. Use the
stage status and output validation recorded for the same run.

## Global checkpoint audit

Every scientific direct command and pipeline stage uses the canonical
`run.resume_policy`. The packaged checkpoint directory is
`run_metadata/checkpoints`; its `checkpoint_events.log` records reuse, warning,
automatic restart, downstream invalidation, and explicit overwrite decisions.
Checkpoint manifests record the resolved-configuration digest, tracked input
and output SHA-256 checksums, software identity, scientific completion status,
and upstream checkpoint dependencies.

With the default `run.resume: true`, PostGWAS reuses only a completely validated
boundary. A real partial run keeps validated earlier stages and reruns the
incomplete stage from its earliest safe boundary. Changed parameters or tracked
inputs generate a visible warning and an automatic restart. Cleanup is allowed
only for checksum-matching PostGWAS-owned regular files inside the configured
output root; modified, unknown, linked, external, or otherwise unproven files
are preserved and stop the restart.

## Failure records

When a stage fails, use the canonical log to identify the failed stage, cause,
corrective guidance, and external-tool output. A partial checkpoint is isolated
from completed results and is never treated as proof of success.

## Related pages

- [Configuration](configuration.md)
- [Pipeline Workflow](pipeline-workflow.md)
- [Input and Output Contracts](input-output-contracts.md)
