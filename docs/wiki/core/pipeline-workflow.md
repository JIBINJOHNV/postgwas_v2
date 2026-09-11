# Pipeline workflow

PostGWAS supports both standalone module commands and a pipeline command. Use a
standalone command when you already have that module's required input artifacts.
Use the pipeline when PostGWAS should run the required preceding steps for one
or more final analyses.

## Harmonisation comes first

Raw summary statistics are prepared with the standalone harmonisation command:

```console
postgwas harmonisation --help
```

Harmonisation is not currently a pipeline-selectable step. Complete it first,
validate the resulting GWAS-VCF, and then provide that VCF to downstream
standalone commands or to `postgwas pipeline`.

## Select final analyses

Pipeline targets describe the results you want, not every prerequisite you
expect to run. For example:

```console
postgwas pipeline --modules finemap --help
```

PostGWAS validates the selected targets and shows their execution order. Shared
preceding steps normally run once. Formatting before imputation and formatting
the resulting VCF for downstream analysis are separate steps because they use
different data.

The contextual help page shows the actual steps and options for the selected
target. Inspect it before starting a run.

## How the execution order is decided

Each target's registered dependencies are expanded recursively and cycles are
rejected. The resulting steps are then emitted in a fixed order:

1. filtering;
2. formatting and imputation, followed by an internal post-imputation filtering
   step when filtering and imputation are both active;
3. LD-block annotation;
4. formatting for downstream tools;
5. the analysis modules;
6. Manhattan plotting;
7. the GWAS-VCF QC summary.

Step directories are numbered by their position in the plan, so the same module
can appear more than once with different numbers, and a different target set
produces different numbers. For `--modules flames` the plan is
`annot_ldblock`, `formatter`, `ld_clump`, `magma`, `magmacovar`, `pops`,
`finemap`, `flames`.

## Optional workflow stages

The pipeline also accepts workflow switches that add filtering, imputation,
plotting, or heritability estimation without naming them as targets:
`--apply-filter`, `--apply-imputation`, `--apply-manhattan`, and
`--heritability`.

```console
postgwas pipeline \
  --modules magma \
  --apply-imputation \
  --apply-manhattan \
  --heritability \
  --help
```

When filtering and imputation are combined, the planner schedules filtering
before imputation and a distinct post-imputation filtering stage. Always inspect
the displayed plan rather than assuming that a command-line option maps to only
one physical step.

`--apply-filter` cannot currently be combined with a target whose plan also
defines a genome-build option — this includes `finemap`, `ld_clump`,
`qc_summary`, `caldera`, `flames`, `mixer`, `gcta_cojo`, and `gcta_gene`. Run
`postgwas sumstat_filter` as a separate step and pass the filtered VCF to the
pipeline instead.

## Export the matching configuration

Generate a configuration for a target and all required preceding steps:

```console
postgwas config export \
  --pipeline finemap \
  --style full \
  --output finemap_pipeline.yaml
```

Review the exported module order and complete all required resource paths before
running the analysis.

## Common input and resource validation

Before scientific stages run, the pipeline checks the entry GWAS-VCF and runs
the registered preflight for each selected module. It collects independent
module failures and presents the observed file checks in a common style.
The YAML audit is written at `pipeline.validation.report_file`; terminal detail
uses one section per file, with indexes grouped. The shared
reporter combines module-specific file facts into that section and shows only
new findings at later stages. Cross-file comparisons refer to stable file
numbers, while stage progress and scientific results remain separate. The
version-2 audit stores each exact file path once with distinct check observations.
The shared
`logging.file_validation.max_screen_files` defaults to 20 successful file
sections and never hides problems; additional successful files are summarized
on screen and retained in full in the saved record set. The audit is updated before
successful completion is announced
or when execution fails, with checks observed at later stages and the failure.

Existing checks have different coverage: some inspect complete file contents,
while others establish availability or header/index validity only. Files produced
by earlier pipeline stages are checked when their consumers can inspect them.
Read [Pipeline Input Validation](pipeline-input-validation.md) before treating
a startup pass as evidence of scientific compatibility.

## Execute the plan

A typical configured downstream run has the following shape:

```console
postgwas pipeline \
  --modules finemap \
  --vcf study.vcf.gz \
  --dataset-id STUDY \
  --output-directory results \
  --run-config finemap_pipeline.yaml
```

The pipeline stops at the failed step if validation, a scientific module, or an
external tool fails. A later step is not reported as complete unless execution
reaches it successfully.

## Related pages

- [Configuration](configuration.md)
- [Input and Output Contracts](input-output-contracts.md)
- [Pipeline Input Validation](pipeline-input-validation.md)
- [Logging and Reproducibility](logging-and-reproducibility.md)
