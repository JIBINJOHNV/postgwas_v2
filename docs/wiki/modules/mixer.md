# MiXeR

## Purpose

MiXeR models single-trait polygenic architecture; GSA-MiXeR tests gene-set
enrichment. The public interface exposes univariate, GSA, or both—not
cross-trait analysis.

## What the analysis does

PostGWAS validates formatter input and complete chromosome reference patterns,
selects native or Docker execution, runs `fit1`/`test1` for univariate analysis,
or splits statistics and fits baseline/full `plsa` models for GSA. It validates
official results, writes summary YAML/TSV, and generates diagnostics if enabled.

## When to use it

Use univariate MiXeR for polygenicity/discoverability and GSA-MiXeR for
model-based gene-set enrichment with build- and ancestry-matched references.

## Input requirements

Formatter `<dataset>_mixer.sumstats.gz`; per-chromosome BIM and LD patterns
containing `@`; MiXeR software or container. GSA also needs annotation patterns
plus baseline, model, and test GO tables, optionally load-library files.

## Command

```console
postgwas mixer --mixer-input-file PATH [--analysis {univariate,gsa,all}] [options]
```

## Minimal example

```console
postgwas mixer \
  --analysis univariate \
  --mixer-input-file formatted/STUDY_mixer.sumstats.gz \
  --bim-file-pattern 'reference/1000G.EUR.QC.@.bim' \
  --ld-file-pattern 'reference/1000G.EUR.QC.@.ld' \
  --dataset-id STUDY \
  --output-directory results
```

## Full example

```console
postgwas mixer \
  --analysis all \
  --mixer-input-file formatted/STUDY_mixer.sumstats.gz \
  --bim-file-pattern 'reference/1000G.EUR.QC.@.bim' \
  --ld-file-pattern 'reference/1000G.EUR.QC.@.ld' \
  --gsa-annotation-file-pattern 'reference/1000G.EUR.QC.@.annot.gz' \
  --gsa-baseline-go-file baseline.tsv \
  --gsa-model-go-file genes.tsv \
  --gsa-test-go-file gene_sets.tsv \
  --genome-build GRCh37 \
  --dataset-id STUDY \
  --output-directory results \
  --seed 10
```

## Parameters

See the generated [Configuration Defaults](../reference/configuration-defaults.md)
for the current analysis mode, backend, genome build, chromosome selection,
resume controls, consistency checks, and GSA settings.

## Processing steps

Validate input schema/chromosomes and `@` patterns; resolve backend; run selected
official stages; validate output identity, convergence, finite estimates, seed
and log evidence; generate diagnostics; write compact summaries.

## Outputs

Univariate raw fit/test JSON and logs are under `results/raw/univariate/`;
summaries are `results/<dataset>_mixer_summary.yaml/.tsv`; figures are under
`plots/`. GSA raw outputs are under `results/raw/gsa/`; summaries are
`<dataset>_gsa_mixer_summary.yaml` and
`<dataset>_gsa_mixer_top_results.tsv`. Logs are run-ID scoped.

## QC and logs

Review variants by chromosome, missing/unexpected chromosomes, optimization
convergence, mixture-versus-infinitesimal AIC/BIC, finite QQ, uncertainty,
model power, seed consistency, diagnostics, and GSA enrichment SE/evidence.

## Interpretation

Polygenicity/discoverability are model estimates conditional on LD and QC. GSA
enrichment requires uncertainty and multiple-gene-set interpretation.

## Common problems

Incomplete patterns, build mismatch, missing chromosome, seed mismatch,
native/container path differences, low power, non-convergence, or invalid GO
schema.

## Limitations

Multi-trait MiXeR is not available through this command.

## Scientific references

- [Frei et al. 2019, MiXeR](https://doi.org/10.1038/s41467-019-10310-0)
- [GSA-MiXeR 2024](https://doi.org/10.1038/s41588-024-01771-1)
- [Official MiXeR implementation](https://github.com/precimed/mixer)
