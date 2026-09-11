# GCTA Gene Analysis

## Purpose

GCTA gene analysis runs fastBAT for genes, fixed segments, or custom SNP sets,
or gene-based mBAT-combo using matched GWAS and PLINK LD-reference data.

## What the analysis does

PostGWAS validates method-specific formatter input, reference BED/BIM/FAM,
build/population declarations and annotations; optionally converts GMT pathways
to reference-specific SNP sets; runs the selected GCTA command in an isolated
staging area; validates its schema; and writes a normalized TSV.

## When to use it

Use it for gene/region/set association with GCTA after creating the shared
GCTA `.ma` formatter file.

## Input requirements

`<dataset>_gcta.ma` for fastBAT or mBAT-combo; PLINK reference prefix;
explicit genome build and reference population; gene list for gene-based
fastBAT, mBAT-combo, and GMT conversion; method-specific set input when
applicable; GCTA 1.94.1 or newer.

## Command

```console
postgwas gcta_gene --gcta-input-file PATH --method METHOD [options]
```

## Minimal example

```console
postgwas gcta_gene \
  --method fastbat_gene \
  --gcta-input-file formatted/STUDY_gcta.ma \
  --gcta-reference-prefix reference/1000G_EUR \
  --gene-list genes_grch37.txt \
  --genome-build GRCh37 \
  --gcta-reference-population EUR \
  --dataset-id STUDY \
  --output-directory results
```

## Full example

```console
postgwas gcta_gene \
  --method fastbat_set \
  --gcta-input-file formatted/STUDY_gcta.ma \
  --gcta-reference-prefix reference/1000G_EUR \
  --gene-list genes_grch37.txt \
  --gmt pathways.gmt \
  --genome-build GRCh37 \
  --gcta-reference-population EUR \
  --gene-window-kb 50 \
  --gcta-reference-maf-min 0.01 \
  --fastbat-ld-cutoff 0.9 \
  --dataset-id STUDY \
  --output-directory results
```

## Parameters

Methods are `fastbat_gene`, `fastbat_segment`, `fastbat_set`, and `mbat_combo`.
See the generated [Configuration Defaults](../reference/configuration-defaults.md)
for the current gene window, segment size, reference MAF, fastBAT and mBAT
settings, frequency checks, output controls, and resume behavior.
For GMT conversion, the configured minimum applies to **coordinate-reference
coverage by the GMT**: matched unique genes divided by all genes in the
coordinate reference. **GMT gene mappability**, with all unique GMT genes as
the denominator, is reported but is not thresholded. Consequently, a focused
GMT can fail the default coverage policy even when every one of its genes maps;
lower the configured threshold explicitly when that is the intended analysis.
The pathway audit table reports complete, partial, or absent mapping and both
gene-ID-mappability and analyzable-gene fractions for every pathway.
`mhc.policy` (`include`, `exclude_snps`, `exclude_genes`, or `exclude_both`)
and `chromosomes.exclude` control the tested genomic scope. The packaged MHC
region is selected from the declared genome build, and a custom region requires
all three `--mhc-chrom`, `--mhc-start`, and `--mhc-end` options.
The `reporting` section configures how many lowest-p-value results are displayed,
their p-value precision, nominal/FDR/family-wise alpha levels, correction-column
names, and chromosome coverage mappings.
`html_report.page_size` configures browser pagination, while
`html_report.columns` defines the displayed result columns independently for
each GCTA method. These presentation settings are loaded from YAML and do not
filter scientific results.

## Processing steps

Resolve config/tool/version; validate build/population and file schemas; check
GWAS/reference IDs and alleles/frequencies; prepare native or GMT-derived sets;
construct the method-specific GCTA command; run or dry-run; validate required
result columns and P values; normalize output; publish the staged files.

## Outputs

Raw primary outputs follow GCTA contracts:
`raw/<dataset>.gene.fastbat`, `.seg.fastbat`, `.fastbat`, or
`.gene.assoc.mbat`. The normalized primary table is
`results/<dataset>_<method>_results.tsv`. A standalone searchable report with
all validated rows is written to
`reports/<dataset>_<method>_gcta_gene_report.html`. It includes input/reference
coverage, analysis scope, nominal/Bonferroni/BH-FDR findings, scientific
warnings and interpretation, configured parameters, and links to the native and
normalized results. Logs and resolved config are under
`logs/` and `run_metadata/`; GMT conversion writes a versioned `prepared_sets/`
bundle with mappings, manifest, README, and checksums.

## QC and logs

Review requested/matched variants, allele/frequency discordance, gene/reference
chromosome overlap, set conversion coverage, GCTA version/command, tested units,
invalid P values, component P values, and primary result validation. The terminal
and canonical log summarize tested units, reference-resolved variants, the
nominal, BH-FDR, and Bonferroni counts, ranked associations, and any scientifically
material input or set omissions. The normalized TSV appends adjusted p-values
and significance flags for every result row without filtering. Gene and segment
summaries also warn about shared chromosomes with no tested units in the result.
MHC/chromosome exclusions are reported with full-reference, shared-variant, and
coordinate-defined gene counts. GCTA receives excluded SNPs by their exact BIM
identifiers; gene-based methods receive a scoped coordinate file, while the
original `.ma` is never rewritten.
The completion manifest fingerprints the HTML report separately from the GCTA
scientific checkpoint. A validated resume can add the report to a legacy run or
refresh changed HTML pagination/columns without rerunning GCTA. An externally
modified report is rejected rather than silently replaced.

## Interpretation

fastBAT and mBAT-combo use different statistics. fastBAT aggregates variant
association evidence while accounting for reference LD; mBAT-combo combines
signed mBAT and unsigned fastBAT evidence. Gene/set results depend on LD,
window/set definitions, reference MAF, and population; do not compare methods
as if they were identical tests or interpret an association as causal.

## Common problems

Missing build/population declaration, wrong formatter file for method, invalid
gene list, no GWAS/reference overlap, incompatible GMT symbols, empty sets,
frequency discordance, or incomplete staged output.

## Limitations

Reference genotypes must be scientifically appropriate; no ancestry or build is
inferred. Confirm ancestry, genome build, and variant compatibility before use.

## Scientific references

- [Bakshi et al. 2016, fastBAT](https://doi.org/10.1038/srep32894)
- [Li et al. 2023, mBAT-combo](https://doi.org/10.1016/j.ajhg.2022.12.006)
- [Benjamini and Hochberg 1995, false discovery rate](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)
- [Official GCTA documentation](https://yanglab.westlake.edu.cn/software/gcta/)
- [MSigDB collection guidance](https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp)
