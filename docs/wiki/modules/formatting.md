# Formatting

## Purpose

Formatting reads one harmonised GWAS-VCF once and creates validated,
tool-specific inputs for downstream PostGWAS modules.

## What the analysis does

The formatter extracts a typed canonical table with bcftools, selects configured
variant IDs independently for each target, infers quantitative versus binary design from case/control count values,
applies each target's required-field and numeric checks, transforms statistics,
and writes outputs atomically. It can produce MAGMA, GCTA gene, SuSiE, FINEMAP,
PRED-LD, LDSC, and MiXeR inputs in one run.

## When to use it

Use it when a downstream analysis requires a tabular contract rather than
GWAS-VCF. Pipeline planning inserts it before dependent analyses and may repeat
it after imputation.

## Input requirements

A harmonised, biallelic GWAS-VCF; dataset ID; output directory; `bcftools`; and
at least one requested format from CLI or YAML. Required VCF fields depend on
the chosen target.

## Command

```console
postgwas formatter --vcf PATH [--format FORMAT [FORMAT ...]] [options]
```

## Minimal example

```console
postgwas formatter \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --dataset-id STUDY \
  --output-directory results \
  --format magma ldsc
```

## Full example

```console
postgwas formatter \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --dataset-id STUDY \
  --output-directory results \
  --format magma gcta_gene susie finemap pred_ld ldsc mixer \
  --run-config formatting.yaml \
  --resume \
  --bcftools bcftools
```

## Parameters

`--format` accepts `magma`, `gcta_gene`, `susie`, `finemap`, `pred_ld`, `ldsc`,
and `mixer` plus documented aliases. When omitted, formats come from YAML.
`--variant-id-type rsid` extracts an rsID from the VCF `ID` field;
`--variant-id-type unique` constructs the configured chromosome-position-REF-ALT
identifier. `variant_identifiers.target_types` can set different conventions
for different outputs. MAGMA pipeline mode inspects BIM field 2 and sets only
the MAGMA target automatically.
Resume reuses outputs only after provenance, content, configuration, and
freshness validation; `--overwrite` takes precedence. Output schemas,
filenames, chromosomes, minimum representable P value, and MiXeR QC come from
the generated [Configuration Defaults](../reference/configuration-defaults.md).

## Processing steps

The formatter validates selection, optionally validates a completion manifest,
extracts `N_ALT=1` records once, types numeric fields, selects each target's
configured ID convention, resolves study design, then runs each exporter in
canonical order. Rows missing or violating a
target's required fields are excluded for that target and counted.

## Outputs

- MAGMA: `<dataset>_magma_snp_loc.tsv` and `<dataset>_magma_p_values.tsv`.
- GCTA gene: `<dataset>_gcta.ma`, shared by fastBAT and mBAT-combo.
- SuSiE: `<dataset>_susie.tsv`.
- FINEMAP: `<dataset>_finemap.tsv`.
- PRED-LD: `pred_ld/<dataset>_chr<chromosome>_pred_ld_input.tsv` for configured chromosomes.
- LDSC: `<dataset>_ldsc_input.tsv`.
- MiXeR: `<dataset>_mixer.sumstats.gz`.
- Provenance: `logs/<dataset>_formatter.log`,
  `run_metadata/resolved_config.yaml`, and
  `run_metadata/formatter_completion.yaml`.

## QC and logs

For each target, review the selected ID type, identifier exclusions, rows
in/out, exclusions for missing or invalid values,
P values bounded at the configured numeric minimum, written schema, inferred
study design, sample-size mode, and output fingerprints.

## Interpretation

Allele direction is preserved: GWAS-VCF ALT is the effect allele. FINEMAP gets
minor allele frequency; MAGMA and GCTA mBAT receive total sample size; SuSiE
uses effective N; LDSC uses case/control counts for binary traits or N from NCO
for quantitative traits; MiXeR receives effective N.

## Common problems

No format selected, no usable IDs of the selected type, duplicate selected IDs,
missing VCF fields, incomplete per-variant case/control
counts, stale outputs, a changed resolved config, or zero rows satisfying one
target's scientific contract.

## Limitations

Formatting does not run the external analyses and cannot make an inappropriate
reference panel compatible. Study design is inferred from count values rather
than a separate VCF trait declaration.

## Scientific references

See [Scientific References](../reference/scientific-references.md) and the
documentation for the downstream tool receiving each export.
