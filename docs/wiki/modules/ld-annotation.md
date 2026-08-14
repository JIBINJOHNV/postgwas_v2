# LD Annotation

## Purpose

LD annotation adds population-specific, approximately independent LD-block
labels to a harmonised GWAS-VCF.

## What the analysis does

For each requested population, PostGWAS annotates intervals from a compressed
BED file into a new VCF INFO field named `<POP>_LDblock`, then re-indexes the
VCF before adding the next population.

## When to use it

Use it before region-based LD pruning and when downstream logic needs the
Berisa–Pickrell LD-block label carried with each variant.

## Input requirements

A harmonised GWAS-VCF, matching genome build, `bcftools`, `tabix`, and one BED
file per population named `<build>_<population>_ldetect.bed.gz`. Each BED has
CHROM, START, END, and annotation-label columns.

## Command

```console
postgwas annot_ldblock --vcf PATH --ld-region-dir PATH [options]
```

## Minimal example

```console
postgwas annot_ldblock \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --genome-build GRCh37 \
  --ld-region-dir reference/ld_blocks \
  --ld-block-populations EUR \
  --dataset-id STUDY \
  --output-directory results
```

## Full example

```console
postgwas annot_ldblock \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --genome-build GRCh37 \
  --ld-region-dir reference/ld_blocks \
  --ld-block-populations EUR AFR EAS \
  --dataset-id STUDY \
  --output-directory results \
  --threads 4
```

## Parameters

Canonical `modules.ld_annotation` YAML supplies GRCh37 and populations EUR,
AFR, and EAS. The corresponding CLI options override those values. Compute
options are shared, although annotation is performed sequentially by
population.

## Processing steps

Copy the input to the output VCF, validate each population BED, add the INFO
header and interval annotation with bcftools, atomically replace the previous
population stage, and create/update the tabix index.

## Outputs

`<output>/<dataset>_ldblock.vcf.gz` and its `.tbi` index. INFO tags are
`<POP>_LDblock` and contain the fourth BED column.

## QC and logs

Confirm every requested INFO header exists, expected variants have block
labels, unassigned variants are understood, and the VCF index opens. The
command does not write a dedicated PostGWAS log.

## Interpretation

LD blocks are population- and build-specific approximations. A block label is
not a causal locus and should not be transferred between populations.

## Common problems

Wrong filename pattern, mismatched `chr` naming, BED/build mismatch, missing
population file, or unavailable bcftools/tabix.

## Limitations

The current engine copies the input before completing all resource validation,
processes populations sequentially, and does not produce a summary report of
assigned/unassigned variants.

## Scientific references

- [Berisa and Pickrell 2016, approximately independent LD blocks](https://doi.org/10.1093/bioinformatics/btv546)
