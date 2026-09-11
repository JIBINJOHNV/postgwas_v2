# MAGMA

## Purpose

MAGMA maps GWAS variants to genes, tests gene association, and optionally runs
competitive gene-set analysis with multiple-testing correction. Configured
positional MAGMA, eMAGMA, H-MAGMA, nMAGMA, and chromMAGMA mappings can run alone
or together without pooling their distinct hypothesis families.

## What the analysis does

In pipeline mode, the formatter creates the final MAGMA tables in one in-memory
export while assessing BIM identifier overlap, optionally retaining that ID
intersection, and applying configured chromosome and MHC SNP policies. MAGMA
then consumes those exact files without rereading or rewriting them, applies the
gene window, and runs
MAGMA annotation and gene analysis, validates optional gene sets and
runs competitive testing, then creates corrected and annotated result tables.
Configured one-character and `\s+` table separators use pandas' fast C parser;
other valid multi-character separator regular expressions automatically use its
Python parser and the selected parser is logged.
Before annotation and testing it applies the MAGMA-local chromosome and MHC
scope. Source VCFs, BIM files, gene locations, precomputed annotations, and GMT
files are never modified.

The default scope is `mhc.policy: exclude_both` and
`chromosomes.exclude: [Y, MT]`; chromosome X remains included. The available
MHC policies are:

- `include`: retain MHC SNPs and annotation units;
- `exclude_snps`: remove MHC SNPs before SNP-to-unit annotation and gene tests;
- `exclude_genes`: retain MHC SNPs but remove annotation units whose tested
  interval overlaps the configured MHC interval;
- `exclude_both`: apply both exclusions.

Excluded chromosomes are removed at both SNP and annotation-unit levels. Their
genes therefore do not enter the effective competitive gene-set universe.
MAGMA reads its build-specific, one-based inclusive MHC interval from
`resources.genomes.<build>.regions.mhc`; the three `--mhc-*` options form an
all-or-none run override.

## When to use it

Use it after `formatter --format magma` for gene-level association and pathway
testing with build- and population-matched reference genotypes and annotations.

## Input requirements

Formatter SNP-location and P-value tables; PLINK BED/BIM/FAM LD-reference
prefix; MAGMA gene-location file; optional GMT/native gene-set file;
MAGMA 1.10 or newer; dataset ID and output directory.

Those are direct-mode inputs for positional MAGMA. Pipeline mode starts from a
harmonised GWAS-VCF and runs the formatter automatically; it still requires the
LD reference and gene annotations. Functional mappings additionally require the
resources declared by their selected mapping definitions. A filename does not
establish a resource's build, gene identifier system, or release.

## Command

```text
postgwas magma --snp-location-file PATH --p-value-file PATH [options]
postgwas pipeline --modules magma --help
```

## Direct mode

### Positional gene association

```console
postgwas magma \
  --snp-location-file formatted/STUDY_magma_snp_loc.tsv \
  --p-value-file formatted/STUDY_magma_p_values.tsv \
  --magma-ld-reference reference/1000G_EUR \
  --gene-location-file reference/NCBI37.3.gene.loc \
  --dataset-id STUDY \
  --output-directory results
```

### Gene and gene-set association

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
  --magma-memory-per-worker-gb 16 \
  --resolve-variants-to-reference \
  --mhc-policy exclude_both \
  --exclude-chromosomes Y MT \
  --dataset-id STUDY \
  --output-directory results
```

### Functional mappings from prepared definitions

Functional annotation definitions are configuration-only resources, so this
workflow uses the generated YAML alongside CLI study inputs. First follow
[resource preparation](../getting-started/resource-setup.md#prepare-magma-functional-mapping-resources)
with `--output-directory reference/magma/functional_mapping`. Its generated
configuration is `configs/magma_functional_mapping.yaml` below that directory.
The pinned bundle is GRCh37/EUR; these example contexts demonstrate different
methods, not a recommendation to combine unrelated tissues in one study.

```console
postgwas magma \
  --run-config reference/magma/functional_mapping/configs/magma_functional_mapping.yaml \
  --magma-mapping positional emagma_brain_amygdala h_magma_adult_brain n_magma_cortex chrom_magma_ovarian_h3k27ac \
  --primary-magma-mapping positional \
  --snp-location-file formatted/STUDY_magma_snp_loc.tsv \
  --p-value-file formatted/STUDY_magma_p_values.tsv \
  --magma-ld-reference reference/magma/functional_mapping/base/ld_reference/g1000_eur/g1000_eur \
  --gene-location-file reference/magma/functional_mapping/base/gene_locations/NCBI37.3/NCBI37.3.gene.loc \
  --dataset-id STUDY \
  --output-directory results/magma_mappings_direct
```

Use only the definition names appropriate to the prespecified question. To run
one definition, give that name to both `--magma-mapping` and
`--primary-magma-mapping`. The same rule applies in pipeline mode below.

| Example definition | Method and target identifiers | Interpretation |
|---|---|---|
| `emagma_brain_amygdala` | eQTL-informed eMAGMA; Entrez | Gene-association p-values |
| `h_magma_adult_brain` | Chromatin-informed H-MAGMA; Ensembl | Gene-association p-values |
| `n_magma_cortex` | Combined nMAGMA; symbols | Gene-association p-values |
| `chrom_magma_ovarian_h3k27ac` | Regulatory-element chromMAGMA; mixed targets | Minimum-element-p ranking, not calibrated gene significance |

## Pipeline mode

### Positional gene and gene-set association

This positional example uses GRCh37/EUR resources and the packaged NCBI37.3
gene-location declaration. The gene-set identifiers must match a supported
identifier column in the gene-location file; review the reported compatibility
decision rather than assuming that every GMT is compatible.

```console
postgwas pipeline \
  --modules magma \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --magma-ld-reference reference/1000G_EUR \
  --gene-location-file reference/NCBI37.3.gene.loc \
  --gene-set-file reference/pathways.gmt \
  --dataset-id STUDY \
  --output-directory results
```

### Functional mapping analyses

With the same prepared GRCh37/EUR bundle and prespecified mappings:

```console
postgwas pipeline \
  --modules magma \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --run-config reference/magma/functional_mapping/configs/magma_functional_mapping.yaml \
  --magma-mapping positional emagma_brain_amygdala h_magma_adult_brain n_magma_cortex chrom_magma_ovarian_h3k27ac \
  --primary-magma-mapping positional \
  --magma-ld-reference reference/magma/functional_mapping/base/ld_reference/g1000_eur/g1000_eur \
  --gene-location-file reference/magma/functional_mapping/base/gene_locations/NCBI37.3/NCBI37.3.gene.loc \
  --dataset-id STUDY \
  --output-directory results/magma_mappings_pipeline
```

The definitions supply annotation paths, build, population, identifier system,
tissue/context and source provenance; selecting a method name alone does not
supply these resources. The formatter creates the paired MAGMA input tables.
Each mapping retains its own hypothesis family. Downstream consumers requiring
calibrated gene statistics reject a chromMAGMA primary result; selecting it
cannot turn its ranking statistic into a calibrated gene p-value.

## Parameters

See the generated [Configuration Defaults](../reference/configuration-defaults.md)
for the current genome build, population, gene windows, gene model, input
schema, harmonisation controls, correction families, output schemas, and resume
behavior. Export the YAML before creating study-specific overrides.

`--magma-memory-per-worker-gb` controls the scheduler's per-worker memory
budget and overrides `modules.magma.batching.memory_per_process_gb`. The number
of concurrent MAGMA gene-analysis workers is bounded by both `--threads` and
the total `--memory-gb` budget. It does not enforce a process-level memory cap.

### Gene identifiers and source declarations

For a positional reference other than the packaged NCBI37.3/Entrez example,
declare the actual first-column identifier system with
`--magma-positional-gene-id-type` and record its source with
`--magma-positional-source-name`, `--magma-positional-source-version`,
`--magma-positional-source-url`, and `--magma-positional-context`. These options
are available in direct and pipeline mode and override the corresponding
`modules.magma.mapping.definitions.positional` keys. They record metadata; they
do not convert Entrez IDs to Ensembl IDs or make incompatible references match.
`--gene-location-alternate-id-type` describes the optional alternate-ID column,
not the primary identifiers. Build and population must also agree with the
resolved MAGMA configuration and all selected mapping definitions.

## Processing steps

Validate executable/version and every input; consume the exact paired MAGMA
tables returned by formatter in pipeline mode; assess exact BIM-field-2 ID
overlap and fail below the configured threshold; optionally retain only matching
IDs; apply the configured duplicate policy to repeated IDs;
apply configured chromosome and MHC SNP exclusions; write harmonized MAGMA
inputs; create a run-owned scoped annotation and gene-exclusion audit; run
annotation and gene analysis; parse results;
apply configured Bonferroni/FDR corrections; optionally validate gene
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

Pipeline-created shared inputs are retained in the pipeline formatter folder;
direct MAGMA creates validated shared inputs in `01_inputs/`. Each mapping has annotations,
any derived pathway-compatible gene-location reference, batches, and native
MAGMA outputs under
`02_intermediates/<mapping-name>/`; corrected and annotated reports are under
`03_results/<mapping-name>/`. The long provenance catalogue is
`04_comparisons/<dataset>_magma_mapping_comparison.tsv`. Provenance is in
`05_logs/<dataset>_magma.log`; the resolved configuration and validated
completion manifest are in `00_run_metadata/`.
The ordered execution summary shown on screen is written to
`04_reports/<dataset>_magma_pipeline_summary.csv`. The companion
`04_reports/<dataset>_magma_report.html` is the scientific results report: it
summarises variant/reference coverage, exclusions, gene and pathway coverage,
multiple-testing results, complete searchable and sortable result tables,
analysis parameters, provenance, complete-result links, and the workflow audit.
HTML page size and displayed columns come from the canonical `html_report`
configuration.
The final screen summary ranks the configured number of strongest gene and
pathway associations by unadjusted MAGMA p-value and displays the corresponding
configured gene corrections and common pathway corrections. The list lengths come from
`screen_summary.top_gene_rows` and `screen_summary.top_pathway_rows` (10 each
by default). `screen_summary.p_value_significant_digits` controls compact
terminal P-value precision; full-precision values remain available in the
result tables and HTML report. Inferential decisions continue to use the
configured correction methods and threshold.

For mappings backed by a gene-location reference,
`<dataset>_magma_genes_annotated.tsv` preserves every native MAGMA gene-result
column, adds the original reference chromosome, start, end, strand, and
alternate gene identifier, and appends the configured multiple-testing
corrections. MAGMA's `START` and `STOP` describe the tested annotation window;
`GENE_REFERENCE_START` and `GENE_REFERENCE_END` preserve the unexpanded gene
interval from the selected reference.

The shared excluded-variant audit is
`02_intermediates/exclusions/<dataset>_magma_excluded_variants.tsv`. Each
mapping writes its scoped `.genes.annot` plus an excluded-unit table under
`02_intermediates/<mapping-name>/`, and an exact exclusion summary under
`03_results/<mapping-name>/`. Annotated gene-set results report the validated
source membership identifiers, tested members, and untested members; the source
pathway file itself remains unchanged.

The resource preparer writes 88 schema-validated GRCh37/EUR mapping definitions
and an eMAGMA Ensembl-to-Entrez network conversion report under
`emagma/networks_entrez/`.
Its generated YAML is a pipeline-compatible configuration with the definitions
nested under `modules.magma`.

## QC and logs

Review formatter ID type, BIM-field-2 identifier overlap, optional reference-ID
filtering, duplicate resolution, gene-ID overlap, MAGMA commands and
version, MHC policy and interval, per-chromosome/MHC exclusion counts, genes
tested, effective pathway membership, correction families, and staged
publication status.

For positional and nMAGMA gene sets, the configured location records contain
primary gene ID, chromosome, integer start/end, `+`/`-` strand, and an optional
alternate ID. Positional MAGMA keeps compatible column-one identifiers unchanged.
When column one is incompatible but column six passes the threshold, PostGWAS
keeps pathway identifiers unchanged and writes a run-owned gene-location file
with columns one and six swapped. Duplicate alternate IDs follow the configured
`gene_sets.alternate_id_duplicate_policy`, and
`input.alternate_gene_id_type` records their identifier system. When neither
identifier system reaches the configured overlap threshold, the default
`gene_sets.identifier_mismatch_action: skip` warns, omits competitive gene-set
analysis, and preserves gene-association results. Use `error` to stop the whole
MAGMA run instead. Pathway identifiers are never substituted or expanded;
non-positional mappings require identifiers compatible with their supplied
annotation. Annotated results retain both source and effective gene memberships.
The screen, CSV, and HTML use the same ordered 12-stage execution summary.
There is no separate preflight screen block. The summary records the LD
reference, gene reference, original pathway file, quantitative identifier
decision, formatter transformations, BIM overlap, analysis-scope policies,
retained inputs, analyses, corrections, and publication.
It explicitly reports that pathway identifiers remain unchanged and identifies
the selected gene-location column. The validated source pathway file is passed
directly to MAGMA without creating a replacement pathway file.
The formatter tables are not given a redundant cross-file screen check because
both come from the same validated retained frame; each table's schema and row
count are still validated and logged. Column six is PostGWAS metadata, because
the [MAGMA manual](https://ibg.colorado.edu/cdrom2021/Day10-posthuma/magma_session/manual_v1.09a.pdf)
defines four required gene-location columns and an optional fifth strand column.
Compatibility uses the smaller reference/effective gene universe, avoiding a
false failure when a comprehensive GMT contains genes outside a coding reference.

Formatter-produced MAGMA tables are already unique: formatter excludes every
row in a duplicated selected-ID group and records the loss before writing its
paired files. Pipeline MAGMA uses those exact files. For independently prepared
direct-module inputs, `snp_harmonisation.duplicate_policy` is applied after any
optional LD-reference filtering. `lowest_p` (the default) retains the row with
the lowest valid p-value, with original row order breaking an exact p-value tie;
`remove` excludes every row in each duplicated-ID group; and `err` stops before
prepared inputs are written. Conflicting coordinate or allele records always
stop under all three policies. The external MAGMA command still receives
`duplicate=error` as a final safeguard, preventing MAGMA's default automatic
`duplicate=drop` behavior from changing the prepared input.

## Interpretation

Gene P values reflect the selected gene window, model, LD reference, and SNP
mapping. Competitive gene-set results test relative association and require
multiple-testing interpretation. Gene results receive the configured
`multiple_testing.gene_methods` (Bonferroni and BH-FDR by default); gene-set
results receive `multiple_testing.global_methods` (Bonferroni, Šidák, Holm, and
BH-FDR by default). `multiple_testing.primary_method` is `fdr_bh` by default and
is the only correction used to declare gene-set significance on screen. Result
tables retain all other global and named-family adjustments as secondary
sensitivity values and include the configured primary method, primary adjusted
p-value, and primary significance decision on every row. This primary correction
applies within one mapping analysis; it does not account for selection across
several mapping definitions.

`exclude_snps` changes gene statistics because MHC variants are absent before
gene testing. `exclude_genes` removes overlapping tested units from gene output
and from the competitive-test background while leaving non-excluded gene
statistics otherwise unchanged. `exclude_both` is the strict combination.
This explicit policy addresses the MHC's unusually dense genes and complex,
long-range LD without claiming that stock MAGMA automatically removes it; it
does not.

## Common problems

Reference-ID mismatch, wrong gene build, insufficient overlap, malformed
gene-set inputs, missing MAGMA license/binary, old version, or
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
- [Genome Reference Consortium GRCh37 MHC region](https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37)
- [Genome Reference Consortium GRCh38 MHC region](https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38)
- [E-MAGMA](https://doi.org/10.1093/bioinformatics/btab115)
- [H-MAGMA](https://doi.org/10.1038/s41593-020-0603-0)
- [nMAGMA](https://doi.org/10.1093/bib/bbaa298)
- [chromMAGMA](https://doi.org/10.26508/lsa.202201446)
- [Benjamini and Hochberg 1995, false-discovery-rate control](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)
- [Benjamini and Bogomolov 2014, selective inference across hypothesis families](https://doi.org/10.1111/rssb.12028)
