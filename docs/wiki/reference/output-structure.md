# Output Structure

PostGWAS groups outputs by dataset and analysis stage. Exact filenames and
subdirectories are module configuration, so consult the resolved YAML and the
module page instead of reconstructing paths from memory.

## Where files are written

Harmonisation writes one tree per dataset below the output directory:
`<dataset>/harmonisation/` for the merged GWAS-VCFs, logs, rejected variants and
QC evidence, `<dataset>/run_metadata/` for the resolved configuration,
sample-sheet row and executed command, and a top-level `run_metadata/` holding
the run summary and run log.

A pipeline run creates one numbered directory per executed step directly below
`--output-directory`, in the form `NN_<step>`. The number is the step's position
in the plan, so it is not a fixed module identifier: a module that runs twice
receives two numbers, and a different target set renumbers every step. The
directory listing is therefore a record of the plan that actually ran. Some step
names differ from the module name — for example `magmacovar` writes
`NN_magma_covar`, and filtering writes `NN_filter_pre_imp` or
`NN_filter_post_imp` depending on where it sits relative to imputation.

Within a step, modules use a consistent internal layout: `results/` for
normalised tables, `raw/` for native external-tool output, `inputs/` for prepared
tool inputs, `logs/`, and `run_metadata/` for the resolved configuration and any
completion manifest.

## Output classes

- Primary results: harmonised VCFs, filtered data, formatted tables, clumps,
  fine-mapping results, gene results, prioritization scores, enrichment tables,
  cross-trait estimates, or plots.
- Indexes and tool inputs: tabix indexes, converted tables, locus files, LD
  matrices, and external-tool parameter files.
- QC evidence: input and output counts, exclusions by reason, missingness,
  allele-frequency comparisons, validation summaries, and rejected records.
- Provenance: resolved configuration, executed commands, software versions,
  stage metadata, logs, timing, and completion status.
- Failure evidence: actionable error summaries and explicitly incomplete stage
  outputs.

## How to interpret completion

The presence of one result file does not prove a stage completed. Confirm the
canonical log reports success and that required companion files and validations
exist. External programs may create partial files before failing.

## Avoiding accidental overwrite

Use a stable dataset identifier and a deliberate output root. Do not reuse an
output directory for scientifically different inputs or configurations unless
the module explicitly supports validated resume behavior. Keep separate runs
when comparing genome builds, populations, reference releases, or thresholds.

## Archiving a run

Archive the original input manifest, resolved configuration, primary outputs,
QC and rejection reports, canonical logs, software/container identity, and a
resource manifest with checksums. Temporary external-tool files can be omitted
only after confirming they are not required to audit or reproduce the result.

