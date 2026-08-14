# LD Clumping

## Purpose

LD clumping identifies approximately independent association signals and
groups them into genomic loci using population-matched linkage disequilibrium
(LD) information.

## What the analysis does

PostGWAS supports one or both configured methods:

- `region` keeps the smallest-P variant in each population-specific annotated
  LD block.
- `standard` performs two r² clumping passes and merges nearby LD boundaries
  into genomic risk loci.

Every selected method is required to succeed. A failed method makes the command
fail, and analysis outputs created or changed during the attempt are renamed
with a `.partial.<run-id>` suffix.

## When to use it

Use LD clumping after harmonisation to identify approximately independent
signals. The `region` method requires LD-block annotations, normally produced
by [LD Annotation](ld-annotation.md). LD clumping is a required pipeline step
for [Fine Mapping](fine-mapping.md). Use reference data whose ancestry and
genome build match the study.

## Input requirements

Both methods require a harmonised, biallelic GWAS-VCF containing the configured
effect, standard-error, allele-frequency, and LP fields. The genome build in the
VCF header must match `modules.ld_clumping.genome_build`.

The `region` method also requires `%INFO/<POP>_LDblock`. The `standard` method
requires a prepared LD-reference directory containing:

- one tabix-indexed `<population>_chr<chromosome>.ld.gz` file per chromosome
  that contains a variant passing `lead_pvalue`; and
- `ld_reference.yaml`, whose build, population, window, minimum r², column
  order, allele-order flag, and symmetric indexing contract match the run.

The standard reference is deliberately symmetric: every PLINK pair is stored
in both endpoint orders. PLINK 1.9 table reports normally contain each pair
once, while PostGWAS runtime tabix queries index the first endpoint.

## Command

Inspect the current command and accepted options:

```console
postgwas ld_clump --help
```

Export the complete canonical configuration:

```console
postgwas config export \
  --module ld_clumping \
  --style full \
  --output ld_clumping.yaml
```

CLI values override canonical YAML values. Argparse supplies no independent
scientific defaults.

## Minimal example

Run annotated-region pruning only:

```console
postgwas ld_clump \
  --vcf study_ldblock.vcf.gz \
  --clumping-methods region \
  --genome-build GRCh37 \
  --population EUR \
  --dataset-id STUDY \
  --output-directory results
```

## Full example

Run both configured methods:

```console
postgwas ld_clump \
  --vcf study_ldblock.vcf.gz \
  --ld-folder reference/ld \
  --genome-build GRCh37 \
  --population EUR \
  --dataset-id STUDY \
  --output-directory results
```

Run standard r² clumping only:

```console
postgwas ld_clump \
  --vcf study.vcf.gz \
  --clumping-methods standard \
  --ld-folder reference/ld \
  --genome-build GRCh37 \
  --population EUR \
  --dataset-id STUDY \
  --output-directory results
```

## Parameters

The packaged configuration follows FUMA-style definitions: independent
significant variants use `lead_pvalue` and `clump_r2`; candidate GWAS variants
are limited by `candidate_pvalue`; lead variants use `lead_r2`; and LD
boundaries within `merge_distance_bp` are merged. Packaged values are P ≤
5×10⁻⁸, candidate P ≤ 0.05, r² 0.6, lead r² 0.1, and 250 kb. All are
schema-validated and configurable.

When `remove_mhc: true`, the configured Genome Reference Consortium coordinates
are applied and the excluded-row count is logged:

- GRCh37: chromosome 6, 28,477,797–33,448,354.
- GRCh38: chromosome 6, 28,510,120–33,480,577.

`missing_index_action: error` is the safe default. `self_only` is an explicit
opt-in assumption that treats a missing index-SNP LD row as having no partners;
the assumption and consequence are logged.

The worker pool respects global `execution.threads` and `execution.memory_gb`.
Its memory estimate is controlled by schema-validated
`modules.ld_clumping.compute.minimum_worker_memory_gb` and
`modules.ld_clumping.compute.input_memory_multiplier` values.

### Allele handling

The GWAS conversion keeps VCF `REF` as the non-effect allele, VCF `ALT` as the
effect allele, and leaves effect size and effect-allele frequency unchanged.
An allele-order-independent ID is used only for matching GWAS and LD records;
it never changes the scientific allele columns. Duplicate canonical variants
with opposite effect-allele orientation cause a hard error. No automatic allele
flip, beta sign change, or frequency transformation is performed.

### Preparing the standard LD reference

PLINK does not run during `postgwas ld_clump`. It runs only when preparing the
LD reference with `tools/resource_preparation/ld_file_preparation.sh`. For each
requested chromosome, the script executes a PLINK 1.9 command using values
supplied explicitly to the script. In this workflow, `--keep-allele-order` is
always included to preserve the original A1/A2 encoding of a PLINK 1 binary
fileset. PLINK 2 uses different LD commands and is not a drop-in replacement.

Complete preparation example:

```console
tools/resource_preparation/ld_file_preparation.sh \
  --bfile /path/to/reference-prefix \
  --output-dir /path/to/ld \
  --population EUR \
  --genome-build GRCh37 \
  --chromosomes 1-22 \
  --plink plink \
  --bgzip bgzip \
  --tabix tabix \
  --window-kb 250 \
  --ld-window-variants 99999 \
  --minimum-r2 0.05 \
  --minimum-maf 0.00001 \
  --threads 8
```

After PLINK finishes, the script duplicates each pair in reverse endpoint
order, sorts it, bgzip-compresses it, builds the tabix index, and writes the
manifest. PostGWAS validates the manifest before analysis.

## Processing steps

1. Validate the selected methods, GWAS-VCF fields, genome build, population,
   configuration, executables, and method-specific references.
2. Extract the configured biallelic GWAS fields without changing allele
   orientation or effect statistics.
3. For `region`, select the smallest-P variant per annotated LD block and write
   both complete and genome-wide-significant pruned results.
4. For `standard`, process chromosomes in a memory-bounded worker pool, query
   the symmetric LD reference, identify significant variants, and apply the two
   configured r² clumping stages.
5. Define association boundaries, merge loci within the configured distance,
   renumber loci, validate outputs, and finalise the canonical log.

## Outputs

Region outputs include the extracted table, all pruned variants, significant
pruned variants, and the region-method log. Standard outputs include the
formatted table, method log, genomic-risk-locus summary, hierarchy,
independent-signal boundaries, and a lead-cluster table when loci exist. Exact
paths come from `modules.ld_clumping.output_layout`.

## QC and logs

The canonical log and resolved-configuration snapshot record selected methods,
thresholds, genome build, MHC exclusions, reference-manifest properties,
allele policy, per-chromosome counts, commands, outputs, and completion or
failure.

Review the input and retained variant counts, genome-wide-significant,
independent, and lead-variant counts, merged-locus count, missing LD queries,
and method completion status. A chromosome with no variant passing
`lead_pvalue` does not require an LD file. Missing or incompatible LD resources
for a significant chromosome are fatal.

## Interpretation

Lead and independent variants represent association structure under the
selected thresholds and reference panel. A lead variant is not necessarily the
causal variant. Compare the locus summaries with the method-specific log before
using the results for fine-mapping or gene prioritisation.

## Common problems

- The `region` method fails when the selected population's LD-block INFO field
  is absent from the VCF.
- The `standard` method fails when `ld_reference.yaml`, a required chromosome
  file, or its tabix index is missing or incompatible.
- Low study/reference variant overlap commonly indicates mismatched build,
  ancestry, identifiers, or allele representation.
- Opposite effect-allele orientations for duplicate canonical variants are
  rejected rather than silently flipped.

## Limitations

Results depend on ancestry matching, genome build, reference coverage,
P/r²/distance thresholds, and the selected method. Region blocks and
FUMA-style standard loci use different definitions and should not be treated
as interchangeable. LD clumping identifies representative signals; it does
not establish causality.

## Scientific references

- [PLINK 1.9 LD reports and window flags](https://www.cog-genomics.org/plink/1.9/ld)
- [PLINK 1.9 `--keep-allele-order`](https://www.cog-genomics.org/plink/1.9/data)
- [FUMA locus and lead-SNP definitions](https://fuma.ctglab.nl/tutorial)
- [Genome Reference Consortium MHC region](https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38)
- [Berisa and Pickrell LD blocks](https://doi.org/10.1093/bioinformatics/btv546)
