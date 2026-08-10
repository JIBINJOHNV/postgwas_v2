# Fine Mapping

## Purpose

Fine mapping estimates which variants within selected association loci are
most compatible with causal effects using SuSiE-RSS or FINEMAP and matched LD.

## What the analysis does

PostGWAS constructs range- or point-based locus windows, applies an LP inclusion
threshold and optional MHC exclusion, aligns formatted summary statistics with
a PLINK reference, computes locus LD, runs the selected model in parallel, and
combines successful locus results into QC and FLAMES-compatible outputs.

## When to use it

Use it after association loci have been defined and after producing the
method-specific formatter table. Fine mapping is meaningful only when LD,
build, ancestry, alleles, sample size, and locus variants are compatible.

## Input requirements

- A locus TSV: `CHROM START END LP` for `range`, or `CHROM POS LP` for `point`.
- `<dataset>_susie.tsv` for SuSiE or `<dataset>_finemap.tsv` for FINEMAP.
- A matching PLINK BED/BIM/FAM reference prefix.
- SuSiE: PLINK plus the required R/susieR environment.
- FINEMAP: PLINK 2, bgenix, LDstore, and FINEMAP executables.

## Command

```console
postgwas finemap --finemap-method {susie,finemap} [options]
```

## Minimal example

```console
postgwas finemap \
  --finemap-method susie \
  --susie-input-file formatted/STUDY_susie.tsv \
  --locus-file loci.tsv \
  --finemap-ld-reference reference/1000G_EUR \
  --dataset-id STUDY \
  --output-directory results
```

## Full example

```console
postgwas finemap \
  --finemap-method finemap \
  --finemap-in-files formatted/STUDY_finemap.tsv \
  --locus-file loci.tsv \
  --locus-type range \
  --window-kb 500 \
  --lp-threshold 7.3 \
  --finemap-ld-reference reference/1000G_EUR \
  --n-causal-snps 5 \
  --prob-cred-set 0.95 \
  --sss \
  --minimum-memory-per-worker-gb 14 \
  --dataset-id STUDY \
  --output-directory results
```

## Parameters

The canonical YAML defaults to SuSiE, range loci, 1,500-kb flanks, LP 7.3, 14 GB per
worker, and skipping chr6:25–35 Mb. SuSiE defaults to 10 components and 180-s
LD/model timeouts. FINEMAP defaults to SSS, 100,000 iterations, five causal
SNPs, 0.95 credible-set probability, and the other values shown by live help.
CLI help displays these YAML values, and only explicitly supplied CLI options
override them.

## Processing steps

Validate tools/resources and versions, prepare/filter locus windows, harmonize
summary and reference alleles, compute and QC LD matrices, schedule loci within
memory limits, fit the model with timeouts, record each locus status, merge
successful outputs, and create minimal FLAMES handoff files.

## Outputs

Both engines use the same categorized layout: published tables, primary
credible sets, and plots are under `results/`; status tables, warnings, and
scientific audits are under `quality_control/`; configuration, software
versions, and run logs are under `run_metadata/`; and the authoritative
post-overlap handoff is under `downstream_inputs/flames/`. Caches, fitted model
objects, prepared locus inputs, primary handoff files, and engine workspaces are
isolated under `intermediate_files/`. Successful SuSiE worker copies are
removed after the published aggregates validate.

## QC and logs

Review attempted/successful/failed loci, variants matched and removed, allele
flips, LD asymmetry/diagonal/eigenvalue checks, model convergence, credible-set
purity/coverage, timeouts, and software versions. The run can complete with a
mixture of successful and failed loci; never infer complete locus coverage from
the top-level status alone.

## Interpretation

Posterior inclusion probabilities and credible sets are conditional on the
chosen locus, variant set, LD matrix, prior/model, and maximum causal effects.
They do not prove biological causality.

## Common problems

Low summary/reference overlap, allele mismatch, wrong build/population, singular
or non-positive LD, insufficient memory, timeout, no locus passing LP, MHC-only
loci, or missing external executables.

## Limitations

Partial locus success requires manual review. Configuration is unified through
the canonical YAML; explicitly supplied CLI values override matching YAML
fields and are recorded in the resolved run configuration.

## Scientific references

- [Wang et al. 2020, SuSiE](https://doi.org/10.1111/rssb.12388)
- [SuSiE-RSS documentation](https://stephenslab.github.io/susieR/reference/susie_rss.html)
- [Benner et al. 2016, FINEMAP](https://doi.org/10.1093/bioinformatics/btw018)
- [FINEMAP documentation](https://christianbenner.com/)
