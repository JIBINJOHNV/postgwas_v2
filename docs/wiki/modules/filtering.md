# Filtering

## Purpose

Filtering materializes a QC-filtered, bgzip-compressed, indexed GWAS-VCF using
one schema-validated filtering configuration. The same resolved configuration
is used by the standalone command and the pipeline.

## What the analysis does

PostGWAS can filter on negative log10 P-value, minor-allele frequency,
imputation quality, study/reference allele-frequency difference, variant type,
palindromic-allele ambiguity, and a build-specific MHC interval. Missing-value
actions are explicit for every nullable field used by an active rule.

`minimum_neglog10_p` is compared with the configured `FORMAT/LP` field. It is
not a raw P-value: genome-wide significance at P ≤ 5×10⁻⁸ corresponds to
approximately `7.30103`.

The configured `genome_build` must be declared in the VCF header using one of
the corresponding `genome_build_header_tokens`. Filename text is never used to
infer a build. When MHC removal is enabled, `mhc.genome_build` must equal the
analysis build, the header must declare contigs, and the configured chromosome
must match a header contig (the equivalent `6`/`chr6` spelling is accepted and
logged). A missing MHC contig is an error because continuing would silently
retain the region.

The default GRCh37 MHC interval is 6:25,000,000–34,000,000 in one-based,
inclusive VCF coordinates. PostGWAS converts its start to zero-based,
half-open BED coordinates when it invokes bcftools. Change the configured
interval and build together when a different scientific definition is needed.
The packaged broad interval follows the extended-MHC convention used in major
GWAS analyses; it remains configurable because published definitions differ.

## When to use it

Use filtering after harmonisation when you need a reproducible GWAS-VCF subset
for downstream analyses with explicit MAF, INFO, significance, variant-type,
strand-ambiguity, reference-frequency, or MHC rules. Do not use it to repair an
incorrect build, allele orientation, or study-sample mapping.

## Input requirements

- A non-empty harmonised `.vcf` or `.vcf.gz`.
- A VCF header declaring the configured build, active FORMAT/INFO tags, and,
  when MHC removal is active, contigs.
- The configured `bcftools`, `tabix`, and `bash` executables.
- A filename-safe dataset ID and an output directory.
- The correct reference-population INFO tag when AF concordance is active.

## Command

The standalone entry point is `postgwas sumstat_filter`. The same filtering
service is also available through a configured PostGWAS pipeline.

## Minimal example

```console
postgwas sumstat_filter \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --dataset-id STUDY \
  --output-directory results \
  --run-config filtering.yaml
```

Standard `--option=value` syntax is supported. Boolean configuration values can
be explicitly enabled or disabled, for example `--remove-mhc` and
`--no-remove-mhc`. CLI values override YAML only when the option is supplied.

## Full example

```console
postgwas sumstat_filter \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --dataset-id STUDY \
  --output-directory results \
  --run-config filtering.yaml \
  --maf-min 0.01 \
  --info-min 0.7 \
  --minimum-neglog10-p 0 \
  --remove-mhc
```

## Parameters

Export the complete canonical configuration before changing less common
policies or field mappings:

```console
postgwas config export --module filtering --style full --output filtering.yaml
```

Important keys include:

- `genome_build`, `genome_build_header_tokens`, and `mhc`
- `maf_min`, `info_min`, `info_max`, and `minimum_neglog10_p`
- `missing_af_action`, `missing_info_action`, and `missing_pvalue_action`
- `reference_population_tag` and `frequency_difference_max`
- `include_indels`, `remove_palindromic`, and `remove_mhc`
- `vcf_fields`, `output_layout`, `sort_output`, and `report_missing_counts`

VCF tags and output path patterns are validated before execution. In
particular, `reference_population_tag` cannot contain shell syntax, and all
configured outputs must remain below the output directory.

## Processing steps

PostGWAS validates configuration, executables, input metadata, build, contigs,
and active VCF tags before constructing expressions. It then performs a
single-stream reason audit, runs a pipefail-protected bcftools pipeline, creates
the tabix index, counts the result, and reconciles mutually exclusive primary
reasons with the observed number removed.

Filter expressions and paths are shell-quoted. The VCF and index are written to
temporary paths and promoted to their configured final names only after
bcftools and tabix succeed. A failed command therefore leaves no zero-byte
final VCF; the canonical log records `STATUS: FAILED` and temporary VCF/index
files are removed.

## Outputs

With the packaged output layout, dataset `STUDY` and build `GRCh37` produce:

- `STUDY_GRCh37_filtered.vcf.gz`
- `STUDY_GRCh37_filtered.vcf.gz.tbi`
- `qc_summary/STUDY_GRCh37_filter_reason_summary.tsv`
- `logs/STUDY_GRCh37_filter_gwas_vcf_bcftools.log`
- `STUDY_GRCh37_mhc_exclude.bed` when MHC removal is active

## QC and logs

The log records the resolved filtering configuration, resolved executable
paths and versions, active expressions, build/field validation, variant counts,
outputs, runtime, warnings, and final status. An unavailable diagnostic count
is reported as unavailable rather than zero.

## Interpretation

The output contains records passing the conjunction of all active include
rules and none of the active exclusion rules. Independent failure counts can
overlap. Primary-reason counts use deterministic first-failure attribution so
they can be reconciled without claiming that a variant failed only one rule.

## Common problems

- A missing configured FORMAT or INFO tag stops preflight before filtering.
- A build mismatch between the VCF header and `genome_build` stops the run.
- An MHC rule requires the configured chromosome to exist in the VCF header.
- Missing values follow their explicit `missing_*_action`; they are never
  silently treated as passing values.

## Limitations

Filtering does not correct upstream allele, build, or metadata errors. It also
does not select one study/sample from a multi-study GWAS-VCF; provide a VCF whose
sample layout already has the intended filtering semantics.

## Scientific references

- [bcftools filtering expressions and target-file behavior](https://samtools.github.io/bcftools/bcftools.html)
- [GWAS-VCF specification](https://github.com/MRCIEU/gwas-vcf-specification)
- [GWAS Catalog summary-statistics format](https://www.ebi.ac.uk/gwas/docs/summary-statistics-format)
- [Sekar et al., extended-MHC analysis at chr6:25–34 Mb](https://www.nature.com/articles/nature16549)
