# Harmonisation Outputs and QC

For dataset `<dataset>`, the canonical dataset root is
`<output>/<dataset>/harmonisation/00_harmonised_sumstat/`. Names below are the
current packaged output-layout patterns; a changed canonical YAML changes the
resolved paths.

## Primary VCFs

The command promises three non-empty primary files:

- `<dataset>_gwas2vcf_<input-build>_merged.vcf.gz`: merged adapter output in the
  inferred input build before the complete annotation/liftover product set;
- `<dataset>_GRCh37_merged.vcf.gz`;
- `<dataset>_GRCh38_merged.vcf.gz`.

The build keys come from `vcf_processing.target_builds`. Companion indexes and
per-chromosome VCF evidence are retained or cleaned according to the canonical
output and cleanup policies.

## Rejections and duplicate evidence

`rejected/` contains per-chromosome and input-level rejected records. The final
`<dataset>_rejected_variants.tsv.gz` consolidates removed records, while
`qc_summary/<dataset>_reject_reasons.tsv` reports every registered reason by
source and total, including explicit zeros. `<dataset>_duplicates.tsv` records
all duplicate-group members, their classification, conflicting fields,
selection status, and action.

Consistent duplicates are selected using configured completeness/quality/N
ordering and a deterministic input-order tie-break. Conflicting association
records are rejected as a group by default; PostGWAS does not select the
smallest P value.

## QC summaries

- `<dataset>_QC_summary.txt`: serialized per-chromosome and dataset processing
  metrics; partial results may be written after chromosome failure.
- `qc_summary/<dataset>_gwas2vcf_summary.tsv`: adapter-input counts.
- `qc_summary/<dataset>_<build>_qc_assessment.tsv`: raw merged-VCF and final
  virtual-subset metrics.
- `qc_summary/<dataset>_<build>_qc_filter_rules.tsv`: independent rule and
  reason counts evaluated against all raw records.
- `qc_summary/<dataset>_<build>_qc_assessment.json`: complete structured QC.

Rule counts can overlap. Only the combined mask defines total excluded and
virtual QC-passed records. Harmonisation does not write a filtered replacement
VCF from this assessment.

## Frequency QC

`eaf_qc/` contains per-chromosome records with missing or out-of-range study
EAF. The main QC reports include study/reference comparable pairs and
frequency-difference summaries where comparison is possible.

## Optional concordance validation

With `--validate`, `qc_summary/concordance/` receives the summary, value
mismatches, input-only variants, VCF-only variants, VCF duplicate records, and
optionally all matched rows. Comparisons orient effects, EAF, and Z statistics
to VCF ALT and can consider direct, swapped, complement, and complement-swapped
alleles under configured policy. Indels require exact representation.

## Logs and provenance

The output root and each dataset have `run_metadata/` with resolved
configuration and command evidence. Within the analysis directory:

- `<dataset>_run_manifest.json` records status, build, VCFs, QC reports, timing,
  chromosome state, and optional validation;
- `logs/<dataset>_dataset.log` is the dataset log;
- `logs/<dataset>_chr<chromosome>.log` preserves chromosome attempts;
- `logs/<dataset>_combined.log` concatenates the final audit trail;
- `logs/adapters/gwas2vcf/` contains exact adapter commands, output, and exit
  codes.

## What to review before downstream use

Confirm final status, inferred input build, complete chromosome set, primary
VCF/index integrity, rejected fractions and reasons, liftover loss, raw and
virtual QC counts, effective-N distribution, study/reference AF concordance,
and any concordance-validation failures. Keep the raw VCF, reports, manifest,
resolved YAML, and logs together.
