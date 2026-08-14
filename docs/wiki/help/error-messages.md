# Error Messages

PostGWAS errors are intended to stop invalid or ambiguous analysis. Fix the
first upstream failure before acting on later missing-output messages.

## Configuration errors

Examples include unknown keys, invalid types/ranges, missing required paths, or
incompatible alternatives. Export the current schema, correct the YAML or
sample sheet, then rerun `postgwas config validate` or `postgwas --validate` as
appropriate. Do not add a guessed fallback.

## Missing or empty file

Check path resolution, container mounts, permissions, file size, compression,
and required companion indexes. Relative harmonisation sample-sheet paths are
relative to the sheet, not the current terminal directory.

## No variants or genes remain

Inspect the preceding QC/rejection report. Typical causes are build/population
mismatch, allele/ID join failure, missing required statistics, strict filters,
no significant loci, or incompatible gene identifiers. Relax a threshold only
after establishing scientific justification.

## External command failed

Use the logged exact command, version, stdout/stderr, and exit code. Confirm the
same command can access all mounted files. Partially created outputs are not a
successful stage.

## Existing or partial output

Every scientific command uses validated resume. A completed checkpoint is reused
only when its parameters, inputs, software, outputs, and checksums still match. A
real partial checkpoint causes a warning and a restart from the earliest safe
module or pipeline-stage boundary; validated earlier stages remain reusable.
Changed parameters or inputs also warn and restart automatically.

PostGWAS replaces only checksum-matching files recorded as its own outputs below
the configured output root. If an output was modified, is a symlink, is unknown,
or lies outside that root, automatic cleanup stops and preserves it for review.
Use `--overwrite` only when forced replacement is intentional, and never combine
files manually across configurations.

## Partial chromosome or locus success

Harmonisation, imputation, clumping, and fine-mapping can report per-unit
status. Confirm whether the top-level command requires all units. A successful
subset is not complete genome coverage; archive the failure list with results.

## Known command limitations

Some commands have important completion limitations, including FLAMES, direct
PoPS, PRED-LD partial coverage, and provider-tolerant enrichment. Read the
relevant module page before interpreting their final messages.
