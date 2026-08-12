# Harmonisation Outputs and QC

For dataset `<dataset>`, the canonical dataset root is
`<output>/<dataset>/harmonisation/`. Names below are the current packaged
output-layout patterns; a changed canonical YAML changes the resolved paths.
This layout applies to new runs; PostGWAS does not relocate results created by
older versions.

## Primary VCFs

The command promises three non-empty primary files:

- `<dataset>_gwas2vcf_<input-build>_merged.vcf.gz`: merged adapter output in the
  inferred input build before the complete annotation/liftover product set;
- `<dataset>_GRCh37_merged.vcf.gz`;
- `<dataset>_GRCh38_merged.vcf.gz`.

The build keys come from `vcf_processing.target_builds`. Companion indexes and
per-chromosome VCF evidence are retained or cleaned according to the canonical
output and cleanup policies.

## Run-level dataset summary

`<output>/run_metadata/harmonisation_run_summary.csv` contains exactly one row
for every dataset selected from the sample sheet, in sample-sheet order. The
file is initialized before the first dataset starts and is atomically updated
after each dataset, so preflight-failed, failed, interrupted, partial, running,
and not-yet-run datasets remain visible rather than disappearing from the run
report.

Columns are grouped from left to right in pipeline order: dataset identity and
status, initial input validation, build inference, study-wide decisions,
chromosome harmonisation, population-frequency comparison, final VCF/QC, and
output provenance. Each row includes:

- final dataset status and a single-line failure reason;
- resolved and independently detected genome build, effect type, and P-value
  type, together with the source of each applied decision and explicit booleans
  showing whether effect and P-value types were selected automatically;
- final frequency type, study strand consensus, and closest reference
  population;
- combined strand-orientation counts across completed chromosomes: evaluated,
  matched, forward, forward-swapped, reverse-complement,
  reverse-complement-swapped, reference-unmatched, palindromic-ambiguous, and
  reference-ambiguous variants;
- expected, completed, and summarized chromosome counts, strand-accounting
  status, comparison panel/population, and the exact chromosome reference files
  represented in the totals;
- total input and input-QC-ready counts, including SNP and indel/other totals
  among rows ready for harmonisation; chromosome-harmonised counts; and final
  merged-VCF counts with their separate SNP and indel/other totals;
- QC-passed and QC-failed SNP counts and percentages, together with the number
  of active rules;
- effect-allele-frequency comparable, concordant, mismatched, and missing
  counts, numeric concordance/mismatch percentages, and the applied AF
  difference cutoff;
- the final effect-scale description, pre-VCF invalid-effect removal count, and
  QC-passed missing/invalid Neff, missing imputation-score, and Neff-outlier
  counts; and
- the dataset manifest and persistent screen-report paths.

Input-ready SNPs are rows whose retained effect and other alleles are both one
of A/C/G/T and exactly one base long. Every other retained row is counted as
indel/other, so those two fields reconcile exactly to `input_ready_variants`.
The combined strand fields reuse counters already saved by each completed
chromosome; PostGWAS does not reread the summary-statistics partitions or VCFs
to create this CSV. The additional final-QC fields are loaded from the QC JSON
already produced from the pipeline's single VCF extraction; they do not trigger
another VCF query. Percent columns contain numeric percentage values, not text
with a percent sign. If any expected chromosome is missing, the strand status
is `partial` and its coverage fields show the incomplete denominator. Therefore
a partial total is never presented as a complete whole-study total.

## Screen report

`<dataset>_screen_report.txt` at the dataset harmonisation root is the
persistent plain-text copy of normal harmonisation stdout. It is initialized
before analysis and flushed as output is produced, so messages already emitted
remain available if a later stage fails. In a multi-dataset run, each report
contains the shared run preparation and final summary plus only that dataset's
processing output; it does not duplicate another dataset's chromosome blocks.
The run-level CSV records the exact report path in `screen_report`.

Terminal display is enabled by default. `--hide-screen` suppresses normal
progress on stdout without disabling report creation; `--show-screen`
explicitly restores terminal display. Fatal CLI errors still use stderr, while
preflight failures are also recorded in each selected dataset's report once
the sample sheet has supplied reliable dataset IDs.

## Rejections and duplicate evidence

`rejected/<dataset>_rejected_variants.tsv.gz` is the single row-level rejection
file. It consolidates input-stage and every chromosome-stage rejection. The
consolidated header and exact row count are validated before the source shards
are removed; a failed write or validation leaves every shard intact and fails
dataset finalization.

`qc_summary/<dataset>_reject_reasons.tsv` reports every registered reason by
source and total, including explicit zeros.
`qc_summary/<dataset>_duplicates.tsv` records all duplicate-group members,
their classification, conflicting fields, selection status, and action.

Consistent duplicates are selected using configured completeness/quality/N
ordering and a deterministic input-order tie-break. Conflicting association
records are rejected as a group by default; PostGWAS does not select the
smallest P value.

## QC summaries

- `qc_summary/<dataset>_chromosomewise_harmonisation_metrics.tsv`: flattened
  chromosome-wise and dataset-level processing metrics; partial results may be
  written after chromosome failure.
- `qc_summary/<dataset>_pre_vcf_column_statistics.tsv`: per-chromosome
  statistics for the harmonised columns exported for VCF creation.
- `qc_summary/<dataset>_gwas2vcf_column_mapping.json`: the validated common
  mapping from GWAS-to-VCF field names to exported column positions.
- `qc_summary/<dataset>_<build>_vcf_qc_metrics.tsv`: raw merged-VCF and final
  virtual-subset metrics.
- `qc_summary/<dataset>_<build>_vcf_qc_rule_results.tsv`: each active rule,
  criterion, action, and affected-variant count evaluated against all raw
  records.
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
one-to-one same-position diagnostic pairs, plus optionally all allele-aware
matched rows. Comparisons orient effects, EAF, and Z statistics to VCF ALT,
compare canonical SE directly, and consider direct, swapped, complement, and
complement-swapped alleles under configured policy. Counts and value summaries
are separated into SNP, indel, and other-variant strata. Input-only and VCF-only
counts are report-only; they cause a warning but never fail the value audit.
The common percentage uses the unique allele-aware union, while source-only
percentages use the corresponding source total. Same-position comparison runs
only for allele-unmatched records and is restricted to unambiguous
one-input/one-VCF pairs of the same variant type. It uses orientation-invariant
values and remains diagnostic; type-mismatched and multiallelic positions are
reported without being paired.

Retained palindromic SNPs are included in effect, EAF, and Z concordance when
the run manifest supplies the forward/reverse study consensus that harmonisation
used to orient them. A missing consensus is reported and never replaced by an
allele-frequency strand guess; SE remains directly comparable.

## Logs and provenance

The output root has `run_metadata/` with the run log, resolved configuration,
and multi-dataset CSV described above. Each dataset also has `run_metadata/`
with its normalized sample-sheet row, resolved configuration, command, and
dataset log. Within the harmonisation directory:

- `<dataset>_screen_report.txt` is the readable terminal transcript described
  above;
- `qc_summary/<dataset>_gwas2vcf_input/` archives the compressed
  per-chromosome TSVs passed to GWAS-to-VCF. Their final `strand_action` column
  records each retained variant's allele-orientation action; the positional
  adapter mapping deliberately excludes that audit-only column;
- `<dataset>_run_manifest.json` records status, build, VCFs, QC reports, timing,
  chromosome state, and optional validation;
- `logs/<dataset>_dataset.log` is the dataset log;
- `logs/<dataset>_chr<chromosome>.log` preserves chromosome attempts;
- `logs/<dataset>_combined.log` concatenates the final audit trail;
- `logs/adapters/gwas2vcf/` contains exact adapter commands, output, and exit
  codes.

Each chromosome EAF QC object separates comparable and discordant reference-AF
counts for palindromic and non-palindromic variants. The packaged policy warns
for non-palindromic discordance and rejects palindromic discordance after
consensus orientation. Strand-reference duplicate counts report exact groups,
conflicting-value groups, discarded groups, and removed reference rows.

## What to review before downstream use

Confirm final status, inferred input build, complete chromosome set, primary
VCF/index integrity, rejected fractions and reasons, liftover loss, raw and
virtual QC counts, effective-N distribution, study/reference AF concordance,
and any concordance-validation failures. Keep the raw VCF, reports, manifest,
resolved YAML, and logs together.
