# LDSC Heritability

## Purpose

The `heritability` command estimates single-trait SNP heritability with LDSC on
formatter-created summary statistics, always on the observed scale and
optionally on the liability scale.

## What the analysis does

It runs CBIIT `munge_sumstats.py` with a merge-alleles list and INFO/MAF
thresholds, then calls `ldsc.py --h2` with reference LD scores and regression
weights. If both sample and population prevalence are provided, it makes a
second liability-scale call.

## When to use it

Use it after LDSC formatting to estimate common-variant heritability and assess
the LDSC intercept/ratio. Liability-scale reporting is for appropriate binary
traits with defensible prevalence values.

## Input requirements

- `<dataset>_ldsc_input.tsv` from the formatter.
- An existing HapMap3-style merge-alleles file; it is required by the current
  runner despite the direct help describing it as optional.
- Matching per-chromosome reference LD scores and regression weights.
- Installed `munge_sumstats.py` and `ldsc.py` entry points.

## Command

```console
postgwas heritability --ldsc-input PATH [options]
```

## Minimal example

```console
postgwas heritability \
  --ldsc-input formatted/STUDY_ldsc_input.tsv \
  --merge-alleles reference/w_hm3.snplist \
  --ref-ld-chr reference/eur_w_ld_chr \
  --w-ld-chr reference/eur_w_ld_chr \
  --dataset-id STUDY \
  --output-directory results
```

## Full example

```console
postgwas heritability \
  --ldsc-input formatted/STUDY_ldsc_input.tsv \
  --merge-alleles reference/w_hm3.snplist \
  --ref-ld-chr reference/eur_w_ld_chr \
  --w-ld-chr reference/eur_w_ld_chr \
  --ldsc-minimum-info 0.9 \
  --ldsc-minimum-maf 0.01 \
  --samp-prev 0.20 \
  --pop-prev 0.01 \
  --dataset-id STUDY \
  --output-directory results
```

## Parameters

Runtime defaults are minimum INFO 0.9, minimum MAF 0.01, and no prevalence
values. Both prevalence values are needed to trigger liability-scale analysis.
The reference and weight prefixes are directories to which LDSC appends
chromosome suffixes.

## Processing steps

Resolve all paths and LDSC entry points, munge the formatter table, validate
the `.sumstats.gz`, optionally run liability-scale h², always run observed-scale
h², parse h²/intercept/ratio for terminal display, and retain full tool logs.

## Outputs

With output prefix `<output>/<dataset>`:

- `<dataset>.sumstats.gz` and munge side outputs;
- `<dataset>.ldsc.log`, the PostGWAS combined command transcript;
- `<dataset>_h2.log` and LDSC result side files for observed scale;
- `<dataset>_Liability_scale_h2.log` and side files when prevalence is complete.

## QC and logs

Review SNPs read/merged/retained by munging, mean chi-square, h² estimate and
standard error, intercept, attenuation ratio, and warnings about weak signal.
Check that alleles and sample-size columns match the formatter contract.

## Interpretation

Observed-scale h² is not liability-scale h². Population prevalence is not the
sample case fraction. LDSC estimates can be unstable or uninformative for low-
heritability traits or mismatched LD references.

## Common problems

Missing merge-alleles path, incompatible HapMap3 IDs/alleles, wrong reference
population/build, missing chromosome LD-score files, invalid prevalence, or
LDSC entry points absent from the active environment.

## Limitations

The public command is single-trait heritability only; it does not expose genetic
correlation in this interface. Supply the merge-alleles file even though the
help text describes it as optional, because the command currently requires it.

## Scientific references

- [Bulik-Sullivan et al. 2015, LD Score regression](https://doi.org/10.1038/ng.3211)
- [CBIIT LDSC implementation](https://github.com/CBIIT/ldsc)
