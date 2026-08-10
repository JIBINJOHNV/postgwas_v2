# Logging and reproducibility

A reproducible PostGWAS result must identify the input data, effective
configuration, software commands, reference resources, transformations, and
completion status that produced it. Terminal progress alone is not an audit
record.

## Two output audiences

PostGWAS separates concise terminal progress from detailed file logging:

- terminal output shows ordered stages, important warnings, failures, and the
  final status needed to monitor a run;
- canonical logs retain detailed parameters, decisions, counts, commands,
  outputs, and failure context.

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

## Failure records

When a stage fails, use the canonical log to identify the failed stage, cause,
corrective guidance, and external-tool output. Do not treat partial files as
completed results.

## Related pages

- [Configuration](configuration.md)
- [Pipeline Workflow](pipeline-workflow.md)
- [Input and Output Contracts](input-output-contracts.md)
