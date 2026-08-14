# MAGMA

## Purpose

MAGMA maps GWAS variants to genes, tests gene association, and optionally runs
competitive gene-set analysis with multiple-testing correction. Configured
positional MAGMA, eMAGMA, H-MAGMA, nMAGMA, and chromMAGMA mappings can run alone
or together without pooling their distinct hypothesis families.

## What the analysis does

PostGWAS validates formatter tables and the LD reference, optionally intersects
exact variant IDs with allele-aware coordinate validation, applies the gene window, runs
MAGMA annotation and gene analysis, optionally converts/validates gene sets and
runs competitive testing, then creates corrected and annotated result tables.

## When to use it

Use it after `formatter --format magma` for gene-level association and pathway
testing with build- and population-matched reference genotypes and annotations.

## Input requirements

Formatter SNP-location and P-value tables; PLINK BED/BIM/FAM LD-reference
prefix; MAGMA gene-location file; optional GMT/native gene-set file;
MAGMA 1.10 or newer; dataset ID and output directory.

## Command

```console
postgwas magma --snp-location-file PATH --p-value-file PATH [options]
```

## Minimal example

```console
postgwas magma \
  --snp-location-file formatted/STUDY_magma_snp_loc.tsv \
  --p-value-file formatted/STUDY_magma_p_values.tsv \
  --magma-ld-reference reference/1000G_EUR \
  --gene-location-file reference/NCBI37.3.gene.loc \
  --dataset-id STUDY \
  --output-directory results
```

## Full example

```console
postgwas magma \
  --snp-location-file formatted/STUDY_magma_snp_loc.tsv \
  --p-value-file formatted/STUDY_magma_p_values.tsv \
  --magma-ld-reference reference/1000G_EUR \
  --gene-location-file reference/NCBI37.3.gene.loc \
  --gene-set-file reference/pathways.gmt \
  --window-upstream 35 \
  --window-downstream 10 \
  --gene-model snp-wise=mean \
  --minimum-snp-overlap 0.5 \
  --minimum-gene-id-overlap 0.5 \
  --resolve-variants-to-reference \
  --dataset-id STUDY \
  --output-directory results
```

## Parameters

See the generated [Configuration Defaults](../reference/configuration-defaults.md)
for the current genome build, population, gene windows, gene model, input
schema, harmonisation controls, correction families, output schemas, and resume
behavior. Export the YAML before creating study-specific overrides.

## Processing steps

Validate executable/version and every input; consume the exact paired MAGMA
tables returned by formatter in pipeline mode; optionally retain exact BIM IDs
and fail below the configured overlap threshold; consolidate a repeated,
coordinate/allele-consistent ID by lowest valid P only for independently
prepared direct-module input;
write harmonized MAGMA inputs; run annotation and gene analysis; parse results;
apply configured Bonferroni/FDR corrections; optionally validate/convert gene
sets, run MAGMA with `--set-annot`, and annotate common tested genes; publish staged
outputs only after completion.

Functional mappings either validate a precomputed eQTL/chromatin annotation,
merge positional/Hi-C/eQTL/TOM components for nMAGMA, or test regulatory
elements and apply the published lowest-element-p gene assignment for
chromMAGMA. That minimum element p-value is a ranking statistic, so PostGWAS
does not apply ordinary gene-level p-value correction or report chromMAGMA gene
significance. Exact external-annotation/BIM ID overlap is required.
The installed metadata preserves the upstream target identifiers: eMAGMA uses
Entrez IDs, H-MAGMA uses Ensembl IDs, nMAGMA uses gene symbols, and chromMAGMA
contains mixed target identifiers.

## Outputs

Validated shared inputs are in `01_inputs/`. Each mapping has annotations,
prepared gene sets, batches, and native MAGMA outputs under
`02_intermediates/<mapping-name>/`; corrected and annotated reports are under
`03_results/<mapping-name>/`. The long provenance catalogue is
`04_comparisons/<dataset>_magma_mapping_comparison.tsv`. Provenance is in
`05_logs/<dataset>_magma.log`; the resolved configuration and validated
completion manifest are in `00_run_metadata/`.

The resource preparer writes 88 schema-validated GRCh37/EUR mapping definitions
and an eMAGMA Ensembl-to-Entrez network conversion report under
`emagma/networks_entrez/`.
Its generated YAML is a pipeline-compatible configuration with the definitions
nested under `modules.magma`.

## QC and logs

Review formatter ID type, optional exact BIM-ID matches, coordinate/allele
conflicts, overlap fractions, duplicate resolution, gene-ID overlap, MAGMA commands and
version, genes tested, correction families, and staged publication status.

For positional and nMAGMA gene sets, the configured location records contain
primary gene ID, chromosome, integer start/end, `+`/`-` strand, and an optional
alternate ID. PostGWAS keeps compatible column-one identifiers unchanged; only
when they are incompatible does it translate matching column-six identifiers to
column one, expanding one-to-many mappings and logging all counts. It fails when
neither identifier system reaches the configured overlap threshold. Annotated
results retain both source and effective gene memberships. Column six
is PostGWAS metadata, because the [MAGMA manual](https://ibg.colorado.edu/cdrom2021/Day10-posthuma/magma_session/manual_v1.09a.pdf)
defines four required gene-location columns and an optional fifth strand column.
Compatibility uses the smaller reference/effective gene universe, avoiding a
false failure when a comprehensive GMT contains genes outside a coding reference.

Formatter-produced MAGMA tables are already unique: formatter excludes every
row in a duplicated selected-ID group and records the loss before writing its
paired files. Pipeline MAGMA uses those exact files. MAGMA's own preparation
still rejects conflicting coordinates/alleles and the external command receives
`duplicate=error`. The configured `lowest_p` rule is a safety policy for direct
MAGMA runs supplied with independently prepared, internally consistent duplicate
rows; it does not restore variants excluded by formatter.

## Interpretation

Gene P values reflect the selected gene window, model, LD reference, and SNP
mapping. Competitive gene-set results test relative association and require
multiple-testing interpretation.

## Common problems

Reference ID/allele mismatch, wrong gene build, insufficient overlap,
incompatible GMT identifiers, missing MAGMA license/binary, old version, or
non-empty partial outputs blocking safe resume. MAGMA preflight validation runs
before staging is created, and an old staging tree containing no files is
removed automatically. A staging tree containing any file remains protected
and requires explicit `--overwrite` after review.

## Limitations

The formatter supports configured rsIDs or coordinate-position-REF-ALT unique
IDs. Mixed/custom BIM identifier systems fail pipeline inference. Results
depend on the selected LD reference and gene annotation.

## Scientific references

- [de Leeuw et al. 2015, MAGMA](https://doi.org/10.1371/journal.pcbi.1004219)
- [Official MAGMA documentation](https://cncr.nl/research/magma/)
- [E-MAGMA](https://doi.org/10.1093/bioinformatics/btab115)
- [H-MAGMA](https://doi.org/10.1038/s41593-020-0603-0)
- [nMAGMA](https://doi.org/10.1093/bib/bbaa298)
- [chromMAGMA](https://doi.org/10.26508/lsa.202201446)
