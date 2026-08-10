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

## Optional workflow stages

The pipeline also accepts workflow switches that add filtering, imputation,
plotting, or heritability estimation:

```console
postgwas pipeline \
  --modules finemap \
  --apply-filter \
  --apply-imputation \
  --apply-manhattan \
  --heritability \
  --help
```

When filtering and imputation are combined, the planner can schedule filtering
before imputation and a distinct post-imputation filtering stage. Always inspect
the displayed plan rather than assuming that a command-line option maps to only
one physical step.

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
- [Logging and Reproducibility](logging-and-reproducibility.md)
