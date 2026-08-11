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
complement, and complement-plus-swapped orientations under policy. It reports
each orientation separately and includes them all in the allele-aware matched
total. It orients effect, effect-allele frequency, and Z to VCF ALT, and compares
canonical standard error directly. When input Z is absent, it derives Z from
effect/SE or signed effect and configured P-value scale.

SNPs, indels, and other variants receive separate matching and value summaries.
The common-variant percentage uses the unique allele-aware union across the two
files as its denominator. Input-only and VCF-only records are report-only and
are shown against their respective source totals. Position diagnostics run only
for those allele-unmatched records. Records at the same chromosome and position
are compared only when there is exactly one record on each side and both have
the same variant type. One-to-one type mismatches and multiallelic positions are
explicitly ambiguous and remain unpaired. Diagnostic comparisons use
absolute effect/Z, folded AF and direct SE because allele orientation is unknown;
they cannot establish variant identity and do not fail the audit. Position
totals are unique chromosome-position counts; type-stratified position counts
may overlap when a mixed or multiallelic position contains several record types.

Default harmonisation YAML includes a retained palindromic SNP in effect, EAF,
and Z concordance when the adjacent run manifest records the strong forward or
reverse study consensus that harmonisation used to orient it. The validator
reconstructs aligned/swapped direction from that consensus; it does not infer
strand from the palindromic allele letters or population AF. If the consensus
is unavailable, the orientation-sensitive values remain excluded and the
audit reports that limitation; SE is still comparable. The defaults allow
strand complement, use tight absolute/relative tolerances, permit no matched-
value mismatch or duplicate/invalid VCF record, and report unmatched variants
without treating their number as a failure.

## Outputs

Below `<dataset>/harmonisation/qc_summary/concordance/`:

- `<dataset>_concordance_summary.tsv`;
- `<dataset>_concordance_mismatches.tsv.gz`;
- `<dataset>_input_only.tsv.gz`;
- `<dataset>_vcf_only.tsv.gz`;
- `<dataset>_vcf_duplicate_records.tsv.gz`;
- `<dataset>_same_position_matches.tsv.gz`;
- optional `<dataset>_all_matches.tsv.gz`.

The complete validation log is under the sibling `logs/` directory.

## Interpretation

Validation tests transformation fidelity, not whether the original study itself
was scientifically valid. Input-only and VCF-only records can reflect documented
filtering, multiallelic splitting, or normalization and are therefore reported
rather than failed. Any value mismatch among allele-aware matches requires
investigation before downstream analysis.
