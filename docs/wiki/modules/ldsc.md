# LDSC Heritability

## Purpose

The `heritability` command estimates single-trait SNP heritability with LDSC on
formatter-created summary statistics, always on the observed scale and
optionally on the liability scale.

## What the analysis does

In the pipeline, it first supplies the merge-alleles list to the LDSC formatter
so duplicated rsIDs can be resolved by compatible alleles. It then supplies the
same file to CBIIT `munge_sumstats.py` with INFO/MAF thresholds and calls
`ldsc.py --h2` with reference LD scores and regression weights. It always runs
observed-scale heritability. When population prevalence is provided, it also
runs liability-scale heritability using explicit or formatter-returned sample
prevalence.

## When to use it

Use it after LDSC formatting to estimate common-variant heritability and assess
the LDSC intercept/ratio. Liability-scale reporting is for appropriate binary
traits with defensible prevalence values.

## Input requirements

- `<dataset>_ldsc_input.tsv` from the formatter.
- An existing HapMap3-style merge-alleles file. It is required by both the
  direct `heritability` command and the heritability pipeline.
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
  --info-min 0.9 \
  --maf-min 0.01 \
  --samp-prev 0.20 \
  --pop-prev 0.01 \
  --dataset-id STUDY \
  --output-directory results
```

### Pipeline example

```console
postgwas pipeline \
  --modules heritability \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --merge-alleles reference/w_hm3.snplist \
  --ref-ld-chr reference/eur_w_ld_chr \
  --w-ld-chr reference/eur_w_ld_chr \
  --dataset-id STUDY \
  --output-directory results
```

## Parameters

`src/postgwas/config/defaults/modules/ldsc.yaml` is the only source of LDSC
analysis defaults. Formatter-specific sample-prevalence aggregation remains in
`src/postgwas/config/defaults/modules/formatting.yaml`. The LDSC analysis
configuration matches the repository-pinned CBIIT LDSC 3.0.1 implementation:
minimum INFO 0.9, minimum MAF 0.01, no explicit minimum N, chunk size 5,000,000,
200 jackknife blocks, `.l2.M_5_50` files, and no constrained intercept,
prevalence, explicit two-step cutoff, or explicit maximum chi-square. Covariance,
delete-value, and retained-MAF outputs are off. Explicit CLI options override the
matching YAML values; omitted CLI options do not own separate defaults. Run
`postgwas heritability --help` for the complete override list.
`--intercept-h2 VALUE` constrains the LD-score regression intercept to that
value; it does not set h². Omit it to retain upstream LDSC's default of
estimating the intercept from the data. A constrained intercept cannot be
combined with an explicit `--two-step` estimator cutoff.

`modules.ldsc.population` and `modules.ldsc.genome_build` both default to
`null`. When explicitly supplied through `--ldsc-population` or
`--ldsc-genome-build`, they are validated and recorded as provenance only.
PostGWAS does not infer either declaration from filenames, and neither setting
changes the LDSC command or selects reference files.

Population prevalence activates liability-scale analysis, which also requires
sample prevalence. In a pipeline, when population prevalence is provided and
`--samp-prev` is omitted, the LDSC service uses only the `sample_prev` value
returned by the formatter. It does not read the formatted table or calculate
sample prevalence. An explicit `--samp-prev` overrides the formatter return.
Direct execution must supply both values for liability-scale analysis. If
population prevalence is omitted, PostGWAS runs observed-scale heritability only;
a sample prevalence supplied on its own does not trigger liability conversion.

For binary traits, the formatter first calculates
`N_CASE / (N_CASE + N_CONTROL)` for every valid LDSC output variant. It then
uses `modules.formatting.ldsc_sample_prevalence.aggregation` to reduce those
fractions to one value. The schema permits `median` (default) or `mean`; it
rejects `sum` because variant rows repeatedly describe the study sample rather
than independent participant groups. The formatter records the aggregation,
variant count, range, and final value, and the pipeline passes that returned
value to this LDSC service. Quantitative traits return no sample prevalence.

Reference and weight paths are directories to which LDSC appends chromosome
numbers and the configured suffixes. The pinned upstream chromosome-split reader
requires autosomes 1–22; PostGWAS validates all required files before execution.
`--merge-alleles` is deliberately the same option name in direct formatter,
direct heritability, and pipeline modes. It is optional only for
`postgwas formatter --format ldsc`; analysis modes require it.

## Processing steps

Resolve the effective YAML-plus-CLI configuration once and record it. Validate
all paths, LDSC entry points, and chromosome-split reference files. In pipeline
mode, format after rsID and allele matching against the merge-alleles table;
direct formatter output can
also request this step explicitly. Munge with the identical reference, validate
the `.sumstats.gz`, optionally run liability-scale h², always run observed-scale
h², require parseable h²/intercept results, and retain full tool logs. Commands
run in a module-owned staging directory; result files are published only after
every requested output validates.

## Outputs

With output prefix `<output>/<dataset>`:

- `<dataset>.sumstats.gz` and munge side outputs;
- `logs/<dataset>_ldsc_service.log`, the canonical PostGWAS command and status log;
- `run_metadata/resolved_config.yaml`, the effective configuration;
- `<dataset>_h2.log` and LDSC result side files for observed scale;
- `<dataset>_liability_h2.log` and side files when prevalence is complete.

After successful output validation, PostGWAS displays a concise LDSC summary by
default: observed-scale h², intercept, attenuation ratio, any liability-scale
findings, prevalence provenance, and result/log paths. `--hide-screen` suppresses
the terminal display while preserving the same summary in the configured screen
transcript. The canonical LDSC service log records the validated estimates as
structured `RESULT` entries.

When enabled, `--print-cov` adds `.cov`; `--print-delete-vals` adds `.delete`
and `.part_delete` for each requested scale.

## QC and logs

Review formatter counts for rsIDs absent from the reference, allele mismatches,
and duplicate groups resolved or left ambiguous. Then review SNPs
read/merged/retained by munging, mean chi-square, h² estimate and standard error,
intercept, attenuation ratio, and warnings about weak signal. Check that alleles
and sample-size columns match the formatter contract.

## Interpretation

Observed-scale h² is not liability-scale h². Population prevalence is not the
sample case fraction. LDSC estimates can be unstable or uninformative for low-
heritability traits or mismatched LD references.

## Common problems

Missing merge-alleles path, incompatible HapMap3 IDs/alleles, wrong reference
population/build, missing chromosome LD-score files, invalid prevalence, or
LDSC entry points absent from the active environment. When validated LDSC
outputs already exist, direct mode stops without modifying them and prints an
actionable message: review the listed files, then use `--overwrite` to replace
them or select a different `--output-directory`. Expected user/configuration
errors are reported without a Python traceback; unexpected internal errors keep
their traceback for diagnosis.

## Limitations

The public command is single-trait heritability only; it does not expose genetic
correlation in this interface. A merge-alleles table can distinguish duplicated
rsIDs only when exactly one record has a compatible allele pair; genuine
ambiguity remains an error. Standard LDSC LD-score files do not contain a
machine-readable genome-build or ancestry manifest. PostGWAS records the
declared `modules.ldsc.population` and `modules.ldsc.genome_build` values and
validates file completeness, but the user must verify that LD scores, weights,
merge-alleles list, GWAS build,
and ancestry come from one compatible reference release.

## Scientific references

- [Bulik-Sullivan et al. 2015, LD Score regression](https://doi.org/10.1038/ng.3211)
- [Pinned CBIIT `ldsc.py` defaults](https://github.com/CBIIT/ldsc/blob/6c673952cee74bd5c57aef1555a03b1c015399a0/ldsc.py)
- [Pinned CBIIT `munge_sumstats.py` defaults](https://github.com/CBIIT/ldsc/blob/6c673952cee74bd5c57aef1555a03b1c015399a0/munge_sumstats.py)
