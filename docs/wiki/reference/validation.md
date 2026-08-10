# Validation Reference

`postgwas --validate` independently compares one original summary-statistics
dataset with its same-build, unfiltered merged GWAS-VCF.

## Required inputs

- The original version-2 harmonisation sample sheet.
- Exact, case-sensitive dataset ID selecting one row.
- The corresponding unfiltered merged GWAS-VCF in the inferred input build.
- Optional harmonisation run YAML controlling tolerances and allele policy.

## Command

```console
postgwas --validate \
  --sample-sheet studies.csv \
  --dataset-id STUDY \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --run-config harmonisation.yaml \
  --output-directory results
```

## Comparisons

The validator matches coordinates and alleles in direct, swapped, strand-
complement, and complement-plus-swapped orientations under policy. It orients
effect, effect-allele frequency, and Z to VCF ALT before comparing. When input Z
is absent, it derives Z from effect/SE or signed effect and configured P-value
scale. Indels are compared exactly and are not FASTA-normalized.

Default harmonisation YAML excludes palindromic comparisons, allows strand
complement, uses tight absolute/relative tolerances, permits no value mismatch,
VCF-only records, or VCF duplicate records, and has no minimum retained fraction.
Export and review current values instead of relying on this summary.

## Outputs

Below `<dataset>/harmonisation/00_harmonised_sumstat/qc_summary/concordance/`:

- `<dataset>_concordance_summary.tsv`;
- `<dataset>_concordance_mismatches.tsv.gz`;
- `<dataset>_input_only.tsv.gz`;
- `<dataset>_vcf_only.tsv.gz`;
- `<dataset>_vcf_duplicate_records.tsv.gz`;
- optional `<dataset>_all_matches.tsv.gz`.

The complete validation log is under the sibling `logs/` directory.

## Interpretation

Validation tests transformation fidelity, not whether the original study itself
was scientifically valid. Input-only records can be legitimate documented
rejections; use harmonisation reason reports to reconcile them. Any value or
orientation mismatch requires investigation before downstream analysis.

