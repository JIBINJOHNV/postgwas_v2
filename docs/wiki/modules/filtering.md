# Filtering

## Purpose

Filtering materializes a QC-filtered, indexed GWAS-VCF using declared thresholds
for association evidence, allele frequency, imputation quality, study/reference
frequency difference, variant type, strand ambiguity, and the MHC region.

## What the analysis does

PostGWAS constructs bcftools include/exclude expressions, counts failures by
rule, applies all active rules, bgzip-compresses the retained records, creates a
tabix index, and reconciles the final removed count with mutually exclusive
primary failure reasons where the audit succeeds.

## When to use it

Use filtering after harmonisation and before modules that should analyze a
restricted variant set. Record the exact thresholds because different filters
change the scientific question and downstream denominator.

## Input requirements

- A harmonised GWAS-VCF with the FORMAT/INFO fields used by enabled rules.
- `bcftools`, `tabix`, a unique dataset ID, and an output directory.
- An appropriate reference-population AF tag when AF concordance is active.

## Command

```console
postgwas sumstat_filter --vcf PATH [options]
```

## Minimal example

```console
postgwas sumstat_filter \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --dataset-id STUDY \
  --output-directory results \
  --run-config filtering.yaml
```

## Full example

```console
postgwas sumstat_filter \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --dataset-id STUDY \
  --output-directory results \
  --minimum-neglog10-p 2 \
  --minimum-maf 0.01 \
  --minimum-info 0.7 \
  --missing-info-action remove \
  --reference-af-column EUR \
  --maximum-af-difference 0.2 \
  --remove-palindromic \
  --palindromic-af-lower 0.4 \
  --palindromic-af-upper 0.6 \
  --remove-mhc \
  --mhc-chrom 6 \
  --mhc-start 25000000 \
  --mhc-end 34000000
```

## Parameters

Use the generated [Configuration Defaults](../reference/configuration-defaults.md)
page or export the current module YAML before use:

```console
postgwas config export --module filtering --style full --output filtering.yaml
```

The current help text contains stale MAF/INFO and flag-description defaults;
use the exported or resolved YAML values for the run.

## Processing steps

PostGWAS validates inputs/tools, counts input variants, evaluates missingness
and each condition in deterministic order, writes the failure-reason audit,
runs a pipefail-protected bcftools pipeline, creates the VCF index, counts the
result, and checks reason-count reconciliation.

## Outputs

If the input filename contains `GRCh37` or `GRCh38`, the build is appended to
the dataset prefix. Outputs are `<prefix>_filtered.vcf.gz`, its `.tbi` index,
`<prefix>_filter_reason_summary.tsv`, and
`logs/<prefix>_filter_gwas_vcf_bcftools.log`. An MHC exclusion BED is also
written when that rule is active.

## QC and logs

Review variants before/after, each independent condition count, each primary
attributed removal, missing-value actions, reconciliation status, and the full
bcftools expressions/command. An unavailable count is reported as unavailable,
not silently converted to zero.

## Interpretation

The output contains only records passing the conjunction of active rules.
Independent rule counts overlap; primary reasons are ordered attribution for
reconciliation, not proof that a variant failed only one condition.

## Common problems

Missing VCF tags can remove all records under remove policies. A wrong AF
population or build creates misleading discordance. MHC contig spelling is
matched to the VCF, but its coordinates must still match the build.

## Limitations

The standalone filtering interface does not expose every lower-level policy.
Filtering cannot correct upstream allele or metadata errors. Its current
default-help mismatch must be resolved in source separately.

## Scientific references

- [bcftools filtering expressions](https://samtools.github.io/bcftools/bcftools.html)
- [GWAS Catalog summary-statistics format](https://www.ebi.ac.uk/gwas/docs/summary-statistics-format)
