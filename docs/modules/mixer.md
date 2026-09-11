# Single-trait MiXeR and GSA-MiXeR

PostGWAS exposes only one-trait analyses from the official GSA-MiXeR software:

- `analysis: univariate` runs `fit1` and `test1` to estimate polygenicity,
  discoverability, native-scale SNP heritability, residual inflation, and model
  diagnostics.
- `analysis: gsa` runs `split_sumstats`, `plsa --gsa-base`, and
  `plsa --gsa-full` to estimate gene-set enrichment for one trait.
- `analysis: all` runs both workflows from the same formatted trait.

PostGWAS does not expose `fit2`, `test2`, a second-trait input, or cross-trait
GSA options. Protected multi-trait tokens are rejected if supplied through
`fit_arguments` or `test_arguments`.

## Scientific input contract

Formatter creates a gzipped, tab-separated file with the configured MiXeR
mapping. The canonical mapping writes:

```text
SNP  CHR  BP  A1  A2  N  Z
```

`A1` is the harmonised effect allele and `A2` is the other allele. For binary
traits, `N` is effective sample size `4/(1/Ncase + 1/Ncontrol)`. For
quantitative traits, `N` is total sample size. Allele frequency is not an input
to the official `fit1`, `test1`, `split_sumstats`, or `plsa` interfaces.

The input genome build and every BIM, LD, SNP-annotation, and load-library
resource must match. MiXeR does not perform liftover. The official distributed
GSA-MiXeR reference data use GRCh37 coordinates, so a different build requires
a complete compatible reference set rather than relabelling the data.

## Commands

Create the one-trait input once:

```sh
postgwas formatter \
  --vcf STUDY_GRCh37.vcf.gz \
  --dataset-id STUDY \
  --output-directory formatted \
  --format mixer
```

Run univariate architecture analysis:

```sh
postgwas mixer \
  --analysis univariate \
  --mixer-input-file formatted/STUDY_mixer.sumstats.gz \
  --dataset-id STUDY \
  --output-directory results \
  --bim-file-pattern '/reference/1000G.EUR.QC.@.bim' \
  --ld-file-pattern '/reference/1000G.EUR.QC.@.run4.ld'
```

Run GSA-MiXeR for the same single trait:

```sh
postgwas mixer \
  --analysis gsa \
  --mixer-input-file formatted/STUDY_mixer.sumstats.gz \
  --dataset-id STUDY \
  --output-directory results \
  --bim-file-pattern '/reference/1000G.EUR.QC.@.bim' \
  --ld-file-pattern '/reference/1000G.EUR.QC.@.run4.ld' \
  --gsa-annotation-file-pattern '/reference/1000G.EUR.QC.@.annot.gz' \
  --gsa-baseline-go-file '/reference/gsa-mixer-baseline-annot.tsv' \
  --gsa-model-go-file '/reference/gsa-mixer-gene-annot.tsv' \
  --gsa-test-go-file '/reference/gsa-mixer-genesetLOO-annot.tsv'
```

`--gsa-loadlib-file-pattern '/reference/loadlib.@.bin'` may replace direct LD
loading for GSA. The BIM and SNP-annotation patterns remain required. Select
`--analysis all` to run architecture and GSA in one invocation.

Pipeline mode adds formatter and selects its MiXeR export automatically:

```sh
postgwas pipeline \
  --modules mixer \
  --analysis gsa \
  --vcf STUDY_GRCh37.vcf.gz \
  --dataset-id STUDY \
  --output-directory results \
  --run-config mixer_pipeline.yaml
```

## Configuration and validation

Export reloadable canonical configuration with either command:

```sh
postgwas --config --module mixer --style minimal > mixer.yaml
postgwas --config --pipeline mixer --style minimal > mixer_pipeline.yaml
```

Paths, chromosome selection, scientific options, thresholds, table schemas,
output layout, compute settings, executable/container locations, and reporting
policy are all configuration values. CLI values override the same keys; the CLI
does not own independent defaults.

Before starting an expensive run, PostGWAS:

1. streams the formatter artifact once to validate its configured schema and
   chromosome coverage;
2. validates every selected-chromosome BIM and required LD file;
3. for GSA, validates every SNP-annotation or load-library file plus the
   configured schema of all three GO-format tables; and
4. records the resolved configuration, selected backend, commands, inputs,
   outputs, failures, and completion status in the canonical log.

## Progress interpretation

Pipeline progress and MiXeR progress use separate counters. For example,
`1/2 · 50%` on the outer pipeline means that one of two pipeline modules has
completed; it does not mean that the active MiXeR optimization is 50% complete.
The univariate MiXeR workflow reports four validated stages: input/resource
validation, `fit1`, `test1`, and result validation/reporting.

Within `fit1` and `test1`, PostGWAS follows events written by the pinned native
MiXeR log. Chromosome LD loading has an exact denominator from the configured
chromosome list. Optimization instead displays the number of observed completed
cost-function evaluations with an unknown denominator (`?`), because the
upstream differential-evolution and Nelder-Mead optimizers stop according to
convergence and do not expose a trustworthy total in advance. PostGWAS therefore
does not convert elapsed time or configured repeats into a speculative percent.
The command reaches 100% only after its exit status and declared outputs pass
validation; a failed command remains below completion.

The GSA settings explicitly pass the configured exclusion ranges, MAF and LD
hard-pruning values, gene-window extension, all-genes label, optional maximum
Z, log-likelihood method, standard-error samples, seed, threads, and optional
Adam schedule. The module never substitutes an implicit scientific threshold.

## Output structure

With canonical output settings, `analysis: all` produces:

```text
02_mixer/
├── results/
│   ├── raw/
│   │   ├── univariate/
│   │   │   ├── STUDY.fit.json
│   │   │   ├── STUDY.fit.log
│   │   │   ├── STUDY.test.json
│   │   │   └── STUDY.test.log
│   │   └── gsa/
│   │       ├── input/STUDY.chr@.sumstats.gz
│   │       ├── STUDY.baseline.json
│   │       ├── STUDY.baseline.log
│   │       ├── STUDY.baseline.snps.csv
│   │       ├── STUDY.baseline.weights
│   │       ├── STUDY.full.json
│   │       ├── STUDY.full.log
│   │       └── STUDY.full.go_test_enrich.csv
│   ├── STUDY_mixer_summary.yaml
│   ├── STUDY_mixer_summary.tsv
│   ├── STUDY_gsa_mixer_summary.yaml
│   └── STUDY_gsa_mixer_top_results.tsv
├── plots/
├── logs/<UTC-run-id>/STUDY_mixer.log
└── run_metadata/resolved_config.yaml
```

The fixed suffixes following each upstream `--out PREFIX` are MiXeR protocol
invariants. Directories and prefixes remain configurable.

The complete upstream enrichment table is preserved. The compact GSA TSV
contains at most the configured `top_results`, restricted to gene sets whose
official `loglike_aic` is above `evidence_aic_threshold` and ranked by estimated
enrichment. This is model-support screening, not a hypothesis-test p-value.
PostGWAS does not manufacture p-values, false-discovery rates, or significance
claims from AIC. The semantic-to-upstream result-column mapping is declared in
`gsa.result_columns`.

GSA-MiXeR is computationally intensive. The upstream guide reports typical
runtimes of roughly 6–12 hours with 8 cores and approximately 30 GB RAM;
dataset and reference choices can change these requirements.

## Scientific sources

The command order, input schema, resource contracts, sample-size equation,
allele meanings, and GSA interpretation follow the
[official MiXeR/GSA-MiXeR documentation](https://github.com/precimed/mixer) and
the [official GSA-MiXeR hello-world workflow](https://github.com/precimed/mixer#gsa-mixer-hello-world-example).
The method is described by Frei et al. in
[Nature Genetics (2024)](https://doi.org/10.1038/s41588-024-01771-1).
