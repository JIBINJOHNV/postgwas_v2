# Logging and reproducibility

A reproducible PostGWAS result must identify the input data, effective
configuration, software commands, reference resources, transformations, and
completion status that produced it. Terminal progress alone is not an audit
record.

## Two output audiences

PostGWAS separates concise terminal progress from detailed file logging:

- terminal output shows ordered stages, important warnings, failures, and the
  final status needed to monitor a run;
- the append-only `logging.screen_log_file` transcript preserves printable
  standard-output and standard-error diagnostics beneath the run output directory,
  without terminal colour codes, cursor controls or PostGWAS live-bar redraws;
- canonical logs retain detailed parameters, decisions, counts, commands,
  outputs, and failure context.

Terminal display is on by default (`logging.show_screen: true`). Every
scientific module and pipeline accepts `--hide-screen` to suppress terminal
display. Hidden runs continue writing the screen transcript; hiding does not
replace or disable the structured canonical logs. Set
`logging.show_screen: true` in YAML to enable display when a run configuration
has disabled it.

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

## Shared terminal style

All public direct commands and pipelines use one runtime theme from
`logging.terminal_style` in the canonical application YAML. Active operations
use cyan, successful completion green, warnings yellow, failures red and
decisions magenta. Counts, ordinary values and paths use the terminal's normal
foreground colour. Symbols and words carry the meaning even without colour.

Shared key-value fields use `logging.terminal_label_width`, with one value
column across nested summaries. Narrow terminals reserve space for the value;
wrapped values continue beneath that column. Narrow live bars shorten the
operation description when necessary to keep the count, percentage and elapsed
time visible; full descriptions remain in the durable milestones. Progress
reports retain their actual stage counts; styling never estimates extra progress.

A run configuration can override the global theme without editing module code:

```yaml
logging:
  terminal_label_width: 42
  terminal_style:
    color: auto            # auto or never
    styles:
      analysis: bold cyan
      success: bold green
      warning: bold yellow
      error: bold red
      decision: bold magenta
      text: default
```

Omitted roles inherit the packaged YAML; invalid roles or Rich styles fail
configuration validation. `color: never` disables colour for the entire run.
`NO_COLOR` is also respected. Redirected output and saved transcripts remain
plain. Live cursor controls may still be used on an interactive terminal when
colour is disabled; they are never saved in the screen transcript.

External-tool wording and printable diagnostics are retained, including literal
square brackets. Native colour/cursor escape sequences are removed before
display and logging so they cannot override the theme. Unclassified native
lines stay neutral and keep their native layout: PostGWAS does not guess their
severity or treat them as validated results. The combined stdout/stderr stream
preserves each stream's order, not a guaranteed ordering between concurrent
writers. Separate structured tool logs remain the detailed audit source.

Implementation uses the documented Rich [console and colour controls](https://rich.readthedocs.io/en/stable/console.html),
[literal text styling](https://rich.readthedocs.io/en/stable/text.html) and
[live console rendering](https://rich.readthedocs.io/en/stable/live.html).

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

Software identity includes a content digest of the installed PostGWAS package,
including Python and R code and bundled configuration/resources, excluding
regenerable Python bytecode caches. Editing source without changing the package
version therefore invalidates prior checkpoints. Keep source unchanged during
an active analysis: this digest identifies code at the stage boundary, not an
immutable snapshot of code that another process might edit while it runs.

Input SHA-256 values are cached with the file device, inode, size, modification
time, and change time. The hash is reused only when that complete identity still
matches; otherwise PostGWAS reads and hashes the file again. A file that changes
during execution prevents checkpoint publication. This keeps large immutable
reference panels fast on repeated runs while retaining content-checksum
provenance.

Output-directory settings and output-layout mappings are not input resources:
relative names such as `results` must not cause another run's live logs to be
tracked. Explicitly supplied upstream result files remain tracked inputs even
when they reside inside an output directory.

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
