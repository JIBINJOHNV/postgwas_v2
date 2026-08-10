# LD Clumping

## Purpose

LD clumping identifies approximately independent significant variants, lead
variants, and merged genomic risk loci using a population-specific LD reference.

## What the analysis does

The command attempts two analyses independently. Region-based pruning keeps the
strongest LP variant per annotated `<POP>_LDblock` and writes the genome-wide
significant subset. The standard path converts GWAS-VCF, queries chromosome LD
files, performs two r² clumping stages, defines locus boundaries, and merges
nearby loci.

## When to use it

Use it after LD annotation when region-based output is wanted, and before
fine-mapping or locus-based prioritization. Use a reference ancestry matching
the study.

## Input requirements

A harmonised GWAS-VCF; population-specific `<POP>_LDblock` tags for the region
path; tabix-indexed `<population>_chr*.ld.gz` reference files for the standard
path; bcftools/tabix; dataset ID and output directory.

## Command

```console
postgwas ld_clump --vcf PATH --ld-folder PATH [options]
```

## Minimal example

```console
postgwas ld_clump \
  --vcf STUDY_ldblock.vcf.gz \
  --ld-folder reference/ld \
  --population EUR \
  --dataset-id STUDY \
  --output-directory results
```

## Full example

```console
postgwas ld_clump \
  --vcf STUDY_ldblock.vcf.gz \
  --ld-folder reference/ld \
  --population EUR \
  --lead-p 5e-8 \
  --r2-clump 0.6 \
  --r2-lead 0.1 \
  --merge-dist 250000 \
  --dataset-id STUDY \
  --output-directory results \
  --threads 8
```

## Parameters

Direct defaults are EUR, lead P `5e-8`, independent-significant r² 0.6, lead
r² 0.1, and locus merge distance 250,000 bp. The canonical YAML also records a
250-kb window and MHC policy, but only settings accepted by the standalone
command should be assumed active.

## Processing steps

Extract biallelic GWAS fields; run block-based pruning; independently convert to
the standard table, process chromosomes in a balanced order, query LD pairs,
identify significant variants, clump at the two thresholds, define boundaries,
merge close loci, renumber loci, and write aggregate tables.

## Outputs

Region path: `<dataset>_vcf.tsv.gz`, `<dataset>_LDpruned_<POP>.tsv`,
`<dataset>_LDpruned_<POP>_sig.tsv`, and `<dataset>_ld_clump.log`.

Standard path: `<dataset>_formatted.tsv`, conversion log,
`<dataset>_GenomicRiskLoci_Summary.txt`,
`<dataset>_GenomicRiskLoci_Hierarchy.txt`,
`<dataset>_IndSig_Clusters_Boundaries.txt`, and
`<dataset>_Lead_Clusters_LD_Only.txt` when loci exist.

## QC and logs

Review LD annotation coverage, formatted variant count, chromosomes and LD
files processed, significant/independent/lead counts, locus count, and failures
from both paths. The command succeeds when either path returns a result, so
confirm which outputs actually exist. The terminal reports chromosomes with
findings or warnings as they are processed, then groups chromosomes without
significant variants into compact ranges in the final summary. Complete
per-chromosome diagnostics remain in the detailed conversion/clumping log.

## Interpretation

Results depend on P/r²/distance thresholds and reference LD. Lead and
independent variants are representatives of association structure, not proof
of causality.

## Common problems

No population LD-block tags, incorrect LD filename/population, low variant-ID
overlap, missing tabix index, no genome-wide significant variants, or one path
failing while the other completes.

## Limitations

The two approaches have different definitions and are not interchangeable.
The command tolerates one method failing and does not emit one unified
completion manifest.

## Scientific references

- [FUMA GWAS locus definition](https://doi.org/10.1038/s41467-017-01261-5)
- [Berisa and Pickrell LD blocks](https://doi.org/10.1093/bioinformatics/btv546)
