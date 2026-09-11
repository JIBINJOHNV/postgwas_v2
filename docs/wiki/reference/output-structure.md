# Output Structure

PostGWAS groups outputs by dataset and analysis stage. Exact filenames and
subdirectories are module configuration, so consult the resolved YAML and the
module page instead of reconstructing paths from memory.

## Where files are written

Harmonisation writes one tree per dataset below the output directory:
`<dataset>/harmonisation/` for the merged GWAS-VCFs, logs, rejected variants and
QC evidence, `<dataset>/run_metadata/` for the resolved configuration,
sample-sheet row and executed command, and a top-level `run_metadata/` holding
the run summary, combined HTML report, and run log. The dataset HTML report is
stored in the harmonisation result root. These reports describe persisted
completion and QC evidence, not an independent reanalysis of the variants.

A pipeline run creates one numbered directory per executed step directly below
`--output-directory`, in the form `NN_<step>`. The number is the step's position
in the plan, so it is not a fixed module identifier: a module that runs twice
receives two numbers, and a different target set renumbers every step. The
directory listing is therefore a record of the plan that actually ran. Some step
names differ from the module name — for example `magmacovar` writes
`NN_magma_covar`, and filtering writes `NN_filter_pre_imp` or
`NN_filter_post_imp` depending on where it sits relative to imputation.

Within a step, module-owned canonical YAML separates normalised results,
native external-tool output, prepared inputs, logs, and run metadata. MAGMA
numbers those classes as `00_run_metadata` through `05_logs` and places each
mapping under its intermediate and result classes. Other modules retain their
documented unnumbered class names; always use the resolved configuration rather
than reconstructing a path from the class name.

Direct runs do not receive the pipeline's `NN_<step>` wrapper. Their own module
layout is written below the requested output directory. For example, the
[Quick Start](../getting-started/quick-start.md) places the QC report under
`results/qc/`, and its MAGMA pipeline under `results/magma_pipeline/`; the
equivalent direct recipe uses `results/formatted/` and `results/magma_direct/`.

## Output classes

- Primary results: harmonised VCFs, filtered data, formatted tables, clumps,
  fine-mapping results, gene results, prioritization scores, enrichment tables,
  single-trait heritability/polygenic-architecture estimates, or plots.
- Indexes and tool inputs: tabix indexes, converted tables, locus files, LD
  matrices, and external-tool parameter files.
- QC evidence: input and output counts, exclusions by reason, missingness,
  allele-frequency comparisons, validation summaries, and rejected records.
- Provenance: resolved configuration, executed commands, software versions,
  stage metadata, logs, timing, and completion status.
- Failure evidence: actionable error summaries and explicitly incomplete stage
  outputs.

## Find the result and its evidence

There is no universal HTML report or identical directory tree for every module.
The module guide names the required outputs for the selected analysis:

| Analysis | What to inspect first | Detailed output guide |
|---|---|---|
| Harmonisation | Dataset/run HTML, merged VCFs, rejected records and reason matrix | [Harmonisation outputs](../harmonisation/outputs-and-qc.md) |
| QC | Self-contained HTML plus its reconciled machine-readable assessment | [QC outputs](../modules/qc-summary.md#outputs) |
| MAGMA | Gene/mapping results, native outputs, exclusions, and generated report | [MAGMA outputs](../modules/magma.md#outputs) |
| LDSC heritability | Observed/liability result logs, munged input, and canonical findings | [LDSC outputs](../modules/ldsc.md#outputs) |
| MiXeR | Native fit/test results, YAML/TSV summaries, and enabled QQ/power diagnostics | [MiXeR outputs](../modules/mixer.md#outputs) |
| Manhattan | PNG/PDF, exact plotted-point TSV, and R transcript; this module does not produce QQ plots | [Manhattan outputs](../modules/manhattan.md#outputs) |
| Pathway enrichment | The needed providers' tables/network outputs and failures separately | [Enrichment outputs](../modules/pathway-enrichment.md#outputs) |

The pathway-enrichment command runs a fixed provider workflow and tolerates
individual provider failures. Check the provider-specific evidence rather than
interpreting a final command message as confirmation that every service returned
usable results.

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
