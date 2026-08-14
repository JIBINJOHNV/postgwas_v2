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

The genome build is inferred only from the authoritative harmonisation metadata
line `##genome_build=<build>`. Each input must contain exactly one such line,
and its complete value must match a build in `resources.genomes`. Filename text,
`##reference`, contig `assembly`, and liftover metadata are not used as
substitutes. This prevents inherited source-build metadata from overriding the
build explicitly assigned to the merged VCF.

When MHC removal is enabled, filtering selects the matching entry from
`mhc_regions`, requires the header to declare contigs, and requires the
configured chromosome to match a header contig. The equivalent `6`/`chr6`
spelling is accepted and logged. A missing build-specific region or MHC contig
is an error because continuing would silently retain the region.

The packaged one-based, inclusive Genome Reference Consortium MHC intervals are
GRCh37 `6:28,477,797–33,448,354` and GRCh38
`6:28,510,120–33,480,577`. PostGWAS converts the selected start to zero-based,
half-open BED coordinates when it invokes bcftools. These intervals describe
the GRC MHC region, not a broader extended-MHC convention; edit the appropriate
`mhc_regions` entry if an analysis protocol requires a different definition.

## When to use it

Use filtering after harmonisation when you need a reproducible GWAS-VCF subset
for downstream analyses with explicit MAF, INFO, significance, variant-type,
strand-ambiguity, reference-frequency, or MHC rules. Do not use it to repair an
incorrect build, allele orientation, or study-sample mapping.

## Input requirements

- A non-empty harmonised `.vcf` or `.vcf.gz`.
- A VCF header containing exactly one authoritative `##genome_build=<build>`
  declaration, active FORMAT/INFO tags, and,
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
  --minimum-maf 0.01 \
  --minimum-info 0.7 \
  --minimum-neglog10-p 7.30103 \
  --remove-mhc
```

Command-line option names are not always identical to the configuration keys
they set: `--minimum-maf` sets `maf_min`, `--minimum-info` sets `info_min`, and
`--maximum-info` sets `info_max`. Use `postgwas sumstat_filter --help` for the
current option names and `postgwas config export --module filtering` for the
current configuration keys.

## Parameters

Export the complete canonical configuration before changing less common
policies or field mappings:

```console
postgwas config export --module filtering --style full --output filtering.yaml
```

Important keys include:

- `mhc_regions` (separate GRCh37 and GRCh38 intervals)
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

If validation fails before a build can be inferred, the failure is recorded in
`logs/STUDY_filter_preflight.log` instead of fabricating a build-specific name.

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
- A missing, duplicate, or unsupported `##genome_build` declaration stops the run.
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
- [Genome Reference Consortium GRCh37 MHC interval](https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37)
- [Genome Reference Consortium GRCh38 MHC interval](https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38)
