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

The liability call repeats the same deterministic LD-score regression with the
same munged statistics, reference scores, weights, estimator settings, and
jackknife blocks. In the pinned CBIIT implementation, sample and population
prevalence are passed only to the final summary: LDSC multiplies h² and its
standard error by the liability conversion factor after fitting the regression.
The intercept and attenuation ratio are therefore shared regression findings,
not separate observed- and liability-scale estimates. PostGWAS requires both
values to match between the two logs before publishing any scientific output.

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

The `null` values for `minimum_n`, `two_step`, and `chisq_max` mean that
PostGWAS omits `--n-min`, `--two-step`, and `--chisq-max`; it does not calculate
or substitute those values. The pinned LDSC implementation therefore retains
its conditional defaults:

- During munging, when the input contains per-variant `N`, the minimum sample
  size is the pandas-default 90th percentile of `N` divided by 1.5. This is not
  0.67 times the maximum `N` unless the sample sizes are effectively constant.
- For heritability with one reference LD-score annotation and an unconstrained
  intercept, LDSC uses its two-step estimator with a chi-square cutoff of 30.
- For heritability with multiple reference LD-score annotations, LDSC does not
  use that automatic two-step setting. When no explicit maximum is supplied,
  it uses `max(0.001 × max(N), 80)` as the maximum chi-square statistic.

These rules are implemented by the pinned external LDSC code, not duplicated
inside PostGWAS. PostGWAS tests instead protect the required boundary: the YAML
values remain `null`, and the corresponding command-line flags remain absent.

The canonical `use_m_5_50: true` setting selects common-SNP `.l2.M_5_50`
files. The public CLI exposes only the original LDSC `--not-M-5-50` switch,
whose displayed default is derived from that YAML setting. With the packaged
configuration, the switch is not used (`Default: false`) and `.l2.M_5_50` files
remain selected. Supplying `--not-M-5-50` overrides the resolved setting to
false and selects `.l2.M` files.

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
calculated from GWAS-VCF case/control counts and returned by the formatter. It
does not read the formatted table or recalculate sample prevalence. An explicit
`--samp-prev` overrides the GWAS-VCF case fraction.
When an explicit CLI or YAML sample prevalence is used in a pipeline, PostGWAS
always compares it with the GWAS-VCF case fraction before invoking
LDSC. The single pre-run screen block, canonical log, and returned result record
both values, their absolute difference, the configured warning threshold, the
case-fraction aggregation and range, and the GWAS-VCF case/control count ranges.
The final summary does not repeat this block; it retains the warning completion
status and reports the selected sample prevalence with the main findings. A
difference greater than
`modules.ldsc.sample_prevalence_comparison.warning_absolute_difference` produces
a scientific warning; the packaged value is `0.01`, meaning one percentage
point. A difference equal to `0.01` does not exceed the threshold.

The explicit value remains authoritative and is not silently replaced. The run
continues after a warning and liability-scale h² uses the CLI/YAML value. To use
the GWAS-VCF case fraction instead, stop the pre-run analysis and rerun without
`--samp-prev`, also leaving `modules.ldsc.sample_prevalence` unset. The pipeline
option `--samp-prev-warning-threshold` can override the YAML warning threshold.
Direct execution has no formatter result to compare and therefore uses the
provided value after schema validation.

In the one pre-run terminal comparison block, the threshold breach and comparison
status remain yellow warnings. When the difference exceeds the threshold, the
complete `Important — value selected` field uses the shared white-on-red attention
style. It states both prevalence values, identifies which value LDSC will use,
and explains how to retain or replace that value. Screen logs remain plain text
without terminal escape codes.

Direct execution must supply both values for liability-scale analysis. If
population prevalence is omitted, PostGWAS runs observed-scale heritability only;
a sample prevalence supplied on its own does not trigger liability conversion.

When both prevalence values are available, PostGWAS keeps LDSC responsible for
the liability conversion rather than reimplementing its formula. It validates
that the observed and liability logs contain identical intercept and ratio
values, including constrained intercepts, undefined ratios, and LDSC's negative-
ratio message. Any difference, or a missing ratio without a constrained
intercept, stops the run before staged results are published. The terminal and
canonical log therefore report the intercept and ratio once; the liability
result is labelled as a conversion of observed-scale h².

The pinned [CBIIT LDSC logger](https://github.com/CBIIT/ldsc/blob/6c673952cee74bd5c57aef1555a03b1c015399a0/ldsc.py#L67-L81)
is designed to duplicate each report line to stdout and the configured `.log`.
Some supported executions can nevertheless exit successfully with a complete
stdout report and a zero-byte log because the upstream file handle is not
explicitly flushed or closed. In exactly that condition, PostGWAS saves the
captured stdout as the staged LDSC log, applies the same strict h², standard
error, intercept, and ratio validation, and records the recovery in the
canonical log. It still stops for a non-zero command, empty stdout, unparseable
metrics, or any missing optional output explicitly requested by configuration.

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
Each reference chromosome requires its `.l2.ldscore.gz` file and the selected
`.l2.M_5_50` or `.l2.M` file. Each regression-weight chromosome requires only
its `.l2.ldscore.gz` file; LDSC reads SNP counts from the reference files and
does not read `.M` files from the weights path.

`--merge-alleles` is deliberately the same option name in direct formatter,
direct heritability, and pipeline modes. It is optional only for
`postgwas formatter --format ldsc`; analysis modes require it.

## Processing steps

Resolve the effective YAML-plus-CLI configuration once and record it. Validate
all paths, LDSC entry points, and chromosome-split reference files. In pipeline
mode, format after rsID and allele matching against the merge-alleles table;
direct formatter output can also request this step explicitly. Munge with the
identical reference, validate the `.sumstats.gz`, run observed-scale h², and then
run liability-scale h² when prevalence is complete. Require parseable
h²/intercept results and require the intercept and ratio to remain unchanged
across both runs. Commands run in a module-owned staging directory; result files
are published only after every requested output and scientific invariant
validates.

## Outputs

With output prefix `<output>/<dataset>`:

- `<dataset>.sumstats.gz` and munge side outputs;
- `logs/<dataset>_ldsc_service.log`, the canonical PostGWAS command and status log;
- `run_metadata/resolved_config.yaml`, the effective configuration;
- `<dataset>_h2.log` and LDSC result side files for observed scale;
- `<dataset>_liability_h2.log` and side files when prevalence is complete.

After successful output validation, PostGWAS displays a concise LDSC summary by
default: observed-scale h², intercept, attenuation ratio, any liability-scale
conversion, prevalence provenance, conversion check, and result/log paths.
`--hide-screen` suppresses the terminal display while preserving the same
summary in the configured screen transcript. The canonical LDSC service log
records one intercept and ratio, the observed h², and the converted liability h²
as structured results.

When enabled, `--print-cov` adds `.cov`; `--print-delete-vals` adds `.delete`
and `.part_delete` for each requested run. LDSC writes these regression and
jackknife files before applying prevalence in its summary method. The observed-
and liability-prefixed copies are therefore shared regression-scale artifacts,
not separately converted liability-scale files.

## QC and logs

Review formatter counts for rsIDs absent from the reference, allele mismatches,
and duplicate groups resolved or left ambiguous. Then review SNPs
read/merged/retained by munging, mean chi-square, h² estimate and standard error,
intercept, attenuation ratio, and warnings about weak signal. Check that alleles
and sample-size columns match the formatter contract. A pipeline prevalence-
warning means the explicit CLI/YAML value differs from the case fraction
calculated from the formatted data by more than the configured absolute
threshold. The analysis continues with the explicit value. Verify the study
ascertainment and remove the override when the formatter value is the intended
sample prevalence. Comparisons at or below the threshold remain visible as
informational checks rather than scientific warnings.

## Interpretation

Observed-scale h² is not liability-scale h². Population prevalence is not the
sample case fraction. The intercept and attenuation ratio describe the fitted
regression and are not liability converted. LDSC estimates can be unstable or
uninformative for low-heritability traits or mismatched LD references.

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
- [Pinned CBIIT heritability and reference readers](https://github.com/CBIIT/ldsc/blob/6c673952cee74bd5c57aef1555a03b1c015399a0/ldscore/sumstats.py)
- [Pinned CBIIT liability conversion and h² summary](https://github.com/CBIIT/ldsc/blob/6c673952cee74bd5c57aef1555a03b1c015399a0/ldscore/regressions.py)
