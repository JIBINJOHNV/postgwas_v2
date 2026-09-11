# MiXeR

## Purpose

MiXeR models single-trait polygenic architecture; GSA-MiXeR tests gene-set
enrichment. The public interface exposes univariate, GSA, or both—not
cross-trait analysis.

## What the analysis does

PostGWAS validates formatter input and complete chromosome reference patterns,
selects native or Docker execution, and runs the official replicated
`fit1`/full-data `test1`/combine procedure for univariate analysis. For GSA it
splits statistics and fits baseline/full `plsa` models. It validates official
results, writes summary YAML/TSV, and generates diagnostics if enabled.

## When to use it

Use univariate MiXeR for polygenicity/discoverability and GSA-MiXeR for
model-based gene-set enrichment with build- and ancestry-matched references.

## Input requirements

Formatter `<dataset>_mixer.sumstats.gz`; per-chromosome BIM and LD patterns
containing `@`; MiXeR software or container. Univariate analysis also requires
the matching 20 random MAF/LD-pruned SNP lists, normally addressed as
`...rep@.snps`. GSA instead needs annotation patterns plus baseline, model, and
test GO tables, optionally load-library files.

The pipeline instead starts from an indexed, single-sample PostGWAS-harmonised
VCF and produces the formatter table automatically. All references must match
the study build and population, and each pattern must cover the configured
chromosomes. In BIM, LD, annotation, and load-library patterns, `@` is the
chromosome placeholder; in the fit-extract pattern it is the replicate
placeholder. It is not a shell wildcard.
For the exact `SNP CHR BP A1 A2 N Z` contract, sample-size semantics and resource
compatibility, see the [MiXeR runtime and resource details](https://github.com/JIBINJOHNV/postgwas_v2/blob/main/docs/modules/mixer.md#scientific-input-contract).

| Analysis | Required external references |
|---|---|
| `univariate` | BIM and MiXeR LD patterns; matching random MAF/LD-pruned `rep@.snps` pattern |
| `gsa` | BIM and SNP-annotation patterns; baseline, model and test GO tables; LD pattern or a compatible GSA load-library pattern |
| `all` | Both sets above; LD remains required for univariate fitting even when GSA uses load-library files |

## Direct mode

```text
postgwas mixer --mixer-input-file PATH [--analysis {univariate,gsa,all}] [options]
```

### Univariate architecture

```console
postgwas mixer \
  --analysis univariate \
  --genome-build GRCh37 \
  --mixer-input-file formatted/STUDY_mixer.sumstats.gz \
  --bim-file-pattern 'reference/1000G.EUR.QC.@.bim' \
  --ld-file-pattern 'reference/1000G.EUR.QC.@.ld' \
  --mixer-fit-extract-file-pattern \
    'reference/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep@.snps' \
  --dataset-id STUDY \
  --output-directory results
```

### GSA gene-set analysis

```console
postgwas mixer \
  --analysis gsa \
  --genome-build GRCh37 \
  --mixer-input-file formatted/STUDY_mixer.sumstats.gz \
  --bim-file-pattern 'reference/1000G.EUR.QC.@.bim' \
  --ld-file-pattern 'reference/1000G.EUR.QC.@.ld' \
  --gsa-annotation-file-pattern 'reference/1000G.EUR.QC.@.annot.gz' \
  --gsa-baseline-go-file baseline.tsv \
  --gsa-model-go-file genes.tsv \
  --gsa-test-go-file gene_sets.tsv \
  --dataset-id STUDY \
  --output-directory results/gsa
```

### Architecture and GSA together

```console
postgwas mixer \
  --analysis all \
  --mixer-input-file formatted/STUDY_mixer.sumstats.gz \
  --bim-file-pattern 'reference/1000G.EUR.QC.@.bim' \
  --ld-file-pattern 'reference/1000G.EUR.QC.@.ld' \
  --mixer-fit-extract-file-pattern \
    'reference/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep@.snps' \
  --gsa-annotation-file-pattern 'reference/1000G.EUR.QC.@.annot.gz' \
  --gsa-baseline-go-file baseline.tsv \
  --gsa-model-go-file genes.tsv \
  --gsa-test-go-file gene_sets.tsv \
  --genome-build GRCh37 \
  --dataset-id STUDY \
  --output-directory results \
  --seed 10
```

For GSA-only execution with precomputed load-library files, replace the
`--ld-file-pattern` argument with
`--gsa-loadlib-file-pattern 'reference/loadlib.@.bin'`. Keep BIM, SNP annotation
and all three GO tables. Do not remove the LD pattern from `--analysis all`.

## Pipeline mode

The planner inserts formatter and selects its MiXeR export. Do not pass
`--mixer-input-file` in these VCF-based commands.

### Univariate architecture

```console
postgwas pipeline \
  --modules mixer \
  --analysis univariate \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --genome-build GRCh37 \
  --bim-file-pattern 'reference/1000G.EUR.QC.@.bim' \
  --ld-file-pattern 'reference/1000G.EUR.QC.@.ld' \
  --mixer-fit-extract-file-pattern \
    'reference/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep@.snps' \
  --dataset-id STUDY \
  --output-directory results/mixer_pipeline
```

### GSA gene-set analysis

```console
postgwas pipeline \
  --modules mixer \
  --analysis gsa \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --genome-build GRCh37 \
  --bim-file-pattern 'reference/1000G.EUR.QC.@.bim' \
  --ld-file-pattern 'reference/1000G.EUR.QC.@.ld' \
  --gsa-annotation-file-pattern 'reference/1000G.EUR.QC.@.annot.gz' \
  --gsa-baseline-go-file baseline.tsv \
  --gsa-model-go-file genes.tsv \
  --gsa-test-go-file gene_sets.tsv \
  --dataset-id STUDY \
  --output-directory results/gsa_pipeline
```

### GSA with precomputed load-library files

```console
postgwas pipeline \
  --modules mixer \
  --analysis gsa \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --genome-build GRCh37 \
  --bim-file-pattern 'reference/1000G.EUR.QC.@.bim' \
  --gsa-loadlib-file-pattern 'reference/loadlib.@.bin' \
  --gsa-annotation-file-pattern 'reference/1000G.EUR.QC.@.annot.gz' \
  --gsa-baseline-go-file baseline.tsv \
  --gsa-model-go-file genes.tsv \
  --gsa-test-go-file gene_sets.tsv \
  --dataset-id STUDY \
  --output-directory results/gsa_loadlib_pipeline
```

### Architecture and GSA together

```console
postgwas pipeline \
  --modules mixer \
  --analysis all \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --genome-build GRCh37 \
  --bim-file-pattern 'reference/1000G.EUR.QC.@.bim' \
  --ld-file-pattern 'reference/1000G.EUR.QC.@.ld' \
  --mixer-fit-extract-file-pattern \
    'reference/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep@.snps' \
  --gsa-annotation-file-pattern 'reference/1000G.EUR.QC.@.annot.gz' \
  --gsa-baseline-go-file baseline.tsv \
  --gsa-model-go-file genes.tsv \
  --gsa-test-go-file gene_sets.tsv \
  --dataset-id STUDY \
  --output-directory results/mixer_and_gsa_pipeline
```

## Parameters

```console
postgwas config export --module mixer --style full --output mixer.yaml
```

The same backend selectors apply to every direct and pipeline recipe:

- `--mixer-backend auto` uses the local software when available, otherwise the
  configured container; this is the packaged default.
- `--mixer-backend native` requires a working local installation. Use `--mixer`
  and `--mixer-figures` to identify nonstandard script locations.
- `--mixer-backend docker` uses the configured container and requires a working
  container runtime. Its image, runtime and platform can be overridden through
  the corresponding `--mixer-container-*` options or canonical YAML.

These choices select execution software, not a different statistical analysis.
The run records the resolved backend and software identity. Quote all chromosome
patterns, and ensure container execution can access every input/reference path.

See the generated [Configuration Defaults](../reference/configuration-defaults.md)
for the current analysis mode, backend, genome build, chromosome selection,
resume controls, consistency checks, and GSA settings.

## Processing steps

Validate input schema/chromosomes and `@` patterns; validate every configured
fit SNP list; resolve backend; run selected official stages; validate output
identity, convergence, finite estimates, seed and log evidence; combine
replicates; generate diagnostics; write compact summaries.

The outer pipeline counter reports completed pipeline modules. MiXeR also shows
its own validation, replicated `fit1`, full-data `test1`, combine, and reporting
stages. During `fit1` and `test1`, chromosome LD loading uses the exact
configured chromosome count. The
optimizer shows observed completed cost-function evaluations as `count/?`, not
a percentage, because the upstream convergence-dependent optimizers do not
provide a trustworthy total in advance. The active phase includes the complete
fit sequence announced by MiXeR, such as
`diffevo-fast (1/2: diffevo-fast → neldermead)`. A command reaches 100% only
after its exit status and outputs validate successfully.

## Outputs

Univariate combined fit/test JSON are under `results/raw/univariate/`; raw
replicate JSON/log files are under `results/raw/univariate/replicates/`;
summaries are `results/<dataset>_mixer_summary.yaml/.tsv`; figures are under
`plots/`. GSA raw outputs are under `results/raw/gsa/`; summaries are
`<dataset>_gsa_mixer_summary.yaml` and
`<dataset>_gsa_mixer_top_results.tsv`. Logs are run-ID scoped. See the
[detailed output structure and interpretation](https://github.com/JIBINJOHNV/postgwas_v2/blob/main/docs/modules/mixer.md#output-structure)
for the native files and the distinction between the complete GSA enrichment
table and the evidence-screened compact summary; the latter is not an FDR or
p-value significance table.

## QC and logs

Review variants by chromosome, missing/unexpected chromosomes, optimization
convergence, mixture-versus-infinitesimal AIC/BIC, finite QQ, uncertainty,
model power, seed consistency, diagnostics, and GSA enrichment SE/evidence.

## Interpretation

Polygenicity/discoverability are model estimates conditional on LD and QC. GSA
enrichment requires uncertainty and multiple-gene-set interpretation.

## Common problems

Incomplete patterns, missing replicate SNP lists, applying a reference bundle
from the wrong ancestry/build, missing chromosome, seed mismatch,
native/container path differences, low power, non-convergence, or invalid GO
schema. Twenty replicated fits are deliberately expensive; the upstream v1.3
notes estimate about ten times the aggregate CPU of the old procedure, and two
threads can require more than a day.

## Limitations

Multi-trait MiXeR is not available through this command.

## Scientific references

- [Frei et al. 2019, MiXeR](https://doi.org/10.1038/s41467-019-10310-0)
- [Holland et al. 2020, univariate MiXeR](https://doi.org/10.1371/journal.pgen.1008612)
- [GSA-MiXeR 2024](https://doi.org/10.1038/s41588-024-01771-1)
- [Official MiXeR implementation](https://github.com/precimed/mixer)
- [Official real-data job](https://github.com/precimed/mixer/blob/master/usecases/mixer_real/MIXER_REAL.job)
