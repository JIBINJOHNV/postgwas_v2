# Harmonisation Outputs and QC

For dataset `<dataset>`, the canonical dataset root is
`<output>/<dataset>/harmonisation/`. Names below are the current packaged
output-layout patterns; a changed canonical YAML changes the resolved paths.
This layout applies to new runs; PostGWAS does not relocate results created by
older versions.

## Input completeness and numeric conversion

The terminal and saved screen report show each field separately, with its actual
study column or selected external/fixed source. **Missing in input** counts
original missing tokens across all input rows. **Invalid numeric** counts
non-missing cells changed to null by scalar conversion, using only the
non-missing values assessed after initial missing-field and coordinate filtering.
These denominators can differ; conversion failures are not automatically row
removals. Counts reuse the existing parsers rather than scanning the data again.

An absent Z column is reported as not supplied, with a planned calculation from
harmonised BETA and SE; it never borrows the effect column's missing count.
Missing identifiers use chromosome_position_effect_other for variants retained
at identifier processing. Missing INFO actions reflect the resolved
`info.on_missing` setting. Multi-value INFO parsing and position validation are
identified separately rather than reported as zero scalar conversion failures.
Nonzero percentages below 0.01% display as `<0.01%`.

Successful conversion does not establish valid scientific ranges or statistical
consistency. The detailed field-completeness TSV retains technical lifecycle
information and adds per-column numeric conversion evidence as JSON.

When effect type is inferred automatically, a structured warning shows the study
column, inferred type, usable-value count, zero/negative count and percentage,
and median. The odds-ratio acceptance threshold and median interval are taken
from the resolved configuration, not fixed display defaults. The warning states
the planned transformation and asks the user to verify the study documentation
before declaring `effect_type`. The detailed explanation remains in the canonical
log. This presentation does not change the inference or scientific calculations.

An automatically selected p-value scale has a matching structured warning at
the same dataset-level step. It shows the exact input column, finite numeric
count, excluded missing/non-numeric/non-finite count, observed range, and the
count and percentage above the resolved detection threshold. The displayed
minimum count, minimum percentage and median interval come from the detector's
recorded configuration. The median check is applied only to a `neglog10`
candidate; selecting `raw` does not establish that every cell is a valid
probability. The notice states whether values will remain raw or be converted
from `-log10(P)` using `P = 10^(-input value)`, and asks the user to confirm
`p_value_type` in the sample sheet or `pvalue.type` in YAML. It does not infer
whether the statistical test is one-sided or two-sided.

This warning is shown once when the study applies automatic inference, not
when an explicit sample-sheet/YAML scale is applied or when detection fails.
Excluded values are omitted from inference, not removed at this step. The notice
is returned with the study decisions and saved in the canonical log; the screen
recorder captures the structured block even with terminal display hidden.
Chromosome workers reuse the decision without repeating this notice. Detection,
conversion, missingness and range policies remain unchanged; no additional
data scan is performed. Any tolerated negative values still follow the
configured chromosome range policy, which can reject, null or fail.

## Primary VCFs

Chromosome screen summaries use separate labelled sections in processing order:
effect-scale conversion, strand/allele alignment and EAF, sample size, effect
inputs, p-values, SE, Z, effect QC, INFO, identifiers and final completeness.
The accounting section reports overall removals once. Strand/EAF removals are
explicitly labelled as a combined total; additional AF removals are the drop
after strand alignment. Reference-frequency disagreements kept by a warning
policy are reported as retained observations, not changed frequencies.

The VCF section distinguishes variants rejected by liftover from successfully
lifted variants excluded by the configured swap policy. These counts reuse
validated liftover accounting. If detailed metrics are unavailable, the report
says so instead of assigning an unverified removal cause or assuming zero.
Saved screen reports use the same layout; already-running processes retain the
formatter loaded when they started.

The command promises two non-empty primary files:

- `<dataset>_GRCh37_merged.vcf.gz`;
- `<dataset>_GRCh38_merged.vcf.gz`.

The raw adapter merge
`<dataset>_gwas2vcf_<input-build>_merged.vcf.gz` is still created and validated
because it is a required merge group. On an `OK` run, PostGWAS deletes it and
its `.tbi` or `.csi` index by default. Supply
`--keep_gwas2vcf_intermediate` to retain and return it as a third output. Failed
or partial finalization can leave the intermediate in place because default
cleanup runs only after post-merge status is `OK`. The build keys come from
`vcf_processing.target_builds`; cleanup uses only the exact configured raw path
and cannot match either final VCF or another dataset.

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
- the dataset manifest, persistent screen-report, and standalone HTML-report
  paths.

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

## HTML reports

Every selected sample has a standalone report at
`<output>/<dataset>/harmonisation/<dataset>_harmonisation_report.html`. It
follows the analysis from input validation to final outputs. Use the opening
findings for a quick review, the section navigation to follow a question, and
expandable evidence for the underlying counts and settings. A processing status
of `OK` means the required operation completed; it does not mean every variant
passed every QC rule or that all warnings can be ignored.

The dataset report has eight numbered sections:

| Section | What it explains |
|---|---|
| 1. Results at a glance | Processing outcome, main variant counts, and findings requiring attention. |
| 2. Inputs and initial checks | Study identity, exact input column mappings, selected internal/external/fixed sources, and initial configuration/resource evidence. |
| 3. Whole-study preparation | All eight dataset steps, including input cleaning, build and strand decisions, statistic-type interpretation, and chromosome preparation. |
| 4. Chromosome harmonisation | All sixteen chromosome steps in execution order, with a compact chromosome overview and expandable evidence for each chromosome. |
| 5. Chromosome completion and variant accounting | Completed/failed chromosomes, attempts, rejected and unprocessed rows, reconciliation, and rejection evidence. |
| 6. Merging and final dataset QC | All five post-merge steps, build-specific VCF outputs, population-frequency comparison, cleanup, virtual QC rules, and final validation. |
| 7. Input–VCF concordance | Optional matching and statistical agreement, separated into SNP, indel, and other-variant results when recorded. |
| 8. Final interpretation and downloads | Follow-up actions, output files, effective configuration, resource paths, timing, and reproducibility evidence. |

Each step explains its purpose, inputs, recorded policy, action and outcome.
The opening **Results at a glance** includes prominent cards for the inferred
genome build, detected effect type, detected P-value type, and **Most similar
reference population**, alongside variant-count cards. Detection is reported separately
from the settings actually used in subsection 1.3; a supplied declaration is
not relabelled as a detected result. Unavailable detection remains unrecorded.
The end-to-end summary then provides ten linked subsections:

1. **Study and run:** input identity, declared trait and recorded interpretation,
   version, completion time, duration, chromosome coverage and resource checks.
2. **Input validation:** rows read, each read-stage removal category, retained
   SNPs and retained indels/other variants, with labelled denominators.
3. **Study-wide decisions:** genome build and supporting testable markers,
   inferred versus declared effect/P-value types, SE scale, frequency
   interpretation and the transformations recorded for export.
4. **Strand and allele frequency:** reference panel/population, consensus,
   orientation counts, unmatched and ambiguous variants, study versus external
   EAF, frequency comparisons, and retained versus rejected disagreements.
5. **Sample size and statistics:** count sources and formula, missingness actions,
   recorded sample-size ranges, supplied versus reconstructed BETA/SE/Z,
   p-value clipping, consistency checks, and INFO provenance/corrections.
6. **Export and liftover:** pre-VCF accounting and rejection reasons, adapter and
   normalized counts, source/target builds, liftover failures, successful lifts
   excluded by policy, and merge outcomes.
7. **Reference-frequency QC:** agreement threshold, comparable/missing variants,
   and a table for every recorded reference population showing correlation,
   mean absolute AF difference, comparable variants and missing reference AF.
   Correlations use their recorded eligible comparison set; missingness uses
   all raw merged-VCF records. The selected-reference AF agreement check has
   its own denominator.
8. **Final VCF QC:** assessed build, all-record and SNP-specific pass/fail counts,
   each active rule, and missingness/low-Neff diagnostics.
9. **Input–VCF concordance:** matching and statistical agreement where recorded;
   not-requested, failed, and missing evidence are not presented as passing.
10. **Final takeaway:** unresolved decisions, cautions and links to outputs and
    supporting evidence. Successful processing is not blanket scientific approval.

These summaries reuse saved manifest, chromosome-stage and QC metrics. They do
not scan raw summary statistics or VCFs again, rerun scientific tests, or alter
any variant. Display percentages use the named recorded denominator; tiny
nonzero fractions are not rounded to zero or complete coverage. Each aggregated
chromosome metric retains its own coverage. Sample-size ranges use extrema,
never averaged chromosome medians. Target liftover totals are explicitly
labelled chromosome-stage totals rather than new merged-VCF measurements.

Checks that were unnecessary or not reached remain visible rather than being
silently omitted. Missing evidence is labelled as not recorded, not applicable,
or not requested when that distinction is supported; it is never converted to
a zero count or an invented pass. Per-chromosome details retain the recorded
stage metrics so supplied statistics, derived values, corrections, exclusions
and retained warnings can be distinguished. Coverage is important for partial
runs: totals from completed chromosomes must not be read as whole-study totals
when other chromosomes did not finish.

Genome-build labels distinguish source-build from target-build VCFs. A
successfully lifted variant excluded by a policy is not the same as a variant
that failed liftover. Likewise, a retained reference-frequency disagreement is
not a changed frequency. Strand and EAF removals are not additional independent
totals to sum when one includes the other.

The population result describes allele-frequency similarity among the reference
populations compared, not confirmation of ancestry. Comparisons can have different usable
denominators; consult their recorded eligibility and missing-value evidence.
The virtual-QC section displays recorded rule criteria and results, explains
overlapping failures, and makes clear that the merged VCF is not physically
filtered by this assessment. Low-Neff findings and other diagnostics should be
interpreted according to their configured action rather than assumed to be
automatic exclusions.

`<output>/run_metadata/harmonisation_run_report.html` begins with a
one-row-per-sample overview and provides dataset search plus expandable
complete report sections for every selected sample in sample-sheet order. It
therefore remains useful when some samples are `OK` and others are `PARTIAL`,
`FAILED`, `PREFLIGHT_FAILED`,
`INTERRUPTED`, or not yet run. Each overview row links to the corresponding
standalone report.

Output links are relative to the report when a report location is known, so
moving the output tree together preserves those links. Full recorded paths
remain available for provenance. A missing output is described instead of
offered as a working download; an intermediate is labelled as intentionally
removed only when the recorded cleanup evidence supports that explanation.
External references may still be unavailable on another computer. Keep the
report, manifest, QC tables and final outputs together when sharing results;
recorded full paths can reveal local directory names.

HTML generation is presentation-only. It reuses the same in-memory run-summary
records and already-read dataset manifests used for the CSV, including their
existing chromosome summaries and QC-report paths. It does not query a VCF,
reread summary statistics, or derive new scientific counts. Reports are
written atomically; failure to write a promised HTML report stops reporting
rather than allowing a nominally successful run without it.

## Screen report

`<dataset>_screen_report.txt` at the dataset harmonisation root is the
persistent plain-text copy of normal harmonisation stdout. It is initialized
before analysis and flushed as output is produced, so messages already emitted
remain available if a later stage fails. In a multi-dataset run, each report
contains the shared run preparation and final summary plus only that dataset's
processing output; it does not duplicate another dataset's chromosome blocks.
The run-level CSV records the exact report path in `screen_report`.

Terminal display is enabled by default. `--hide-screen` suppresses the terminal
copy without disabling the dataset reports or the shared transcript. The
packaged transcript path is `<output>/run_metadata/screen.log`, and
`logging.screen_log_file` can change the relative location. Set
`logging.show_screen: true` in YAML to restore display when a run configuration
has disabled it. The shared transcript captures stdout and stderr, including
fatal CLI errors. Preflight failures are also recorded in each selected
dataset's report once the sample sheet has supplied reliable dataset IDs.

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
  scientific category and purpose, criterion, virtual-subset decision,
  input-wide percentage, unique-only impact, overlap, and directional trigger
  counts evaluated against all input records.
- `qc_summary/<dataset>_<build>_qc_assessment.json`: complete structured QC.
- `qc_summary/<dataset>_<build>_qc_summary.csv`: combined overall, metric,
  rule, validation, and execution-provenance rows for portable review.
- `reports/<dataset>_<build>_qc_report.html`: detailed self-contained QC report
  with key findings, PostGWAS AF/score provenance, readable rule cards, inactive
  rules, and exact accounting, rendered from the same reconciled assessment
  without rereading the VCF.

Rule counts can overlap. Only the combined mask defines total excluded and
virtual QC-passed records. Harmonisation does not write a filtered replacement
VCF from this assessment.

The report separates score quality/completeness, reference concordance, allele
compatibility, and analysis-scope decisions. It describes matching variants as
excluded from the virtual subset rather than treating every SNP-only, MHC, or
palindromic policy match as a defective variant. Configured PostGWAS provenance
also states whether the assessed score is study measured or an external proxy.

This assessment is owned by the QC module and uses the schema-validated
`modules.qc_summary` rules, VCF fields, build-specific MHC interval, table
settings, and output layout. Harmonisation invokes that service internally and
includes the resolved QC configuration in its provenance, so `postgwas qc` and
harmonisation do not maintain separate summary policies.

Harmonisation applies this assessment to its input-build merged VCF. The QC
service infers the coordinate build from that VCF's exact
`##genome_build=<build>` declaration and selects the matching output names and
MHC interval; QC does not carry a separate configured target-build value.

Every merged VCF declares its coordinate system explicitly as
`##genome_build=<build>`. Source-build annotated and raw GWAS-to-VCF outputs use
the detected input build; successfully lifted outputs use the target build. A
not-lifted VCF uses the input build because those records never entered the
target coordinate system. PostGWAS adds this line while the existing
chromosome-concatenation stream is written, using uncompressed BCF between
bcftools commands, so no second full-file pass is required. This follows the
official [bcftools streaming guidance](https://samtools.github.io/bcftools/howtos/scaling.html)
and the VCF `##key=value` metadata convention in the
[VCF specification](https://samtools.github.io/hts-specs/VCFv4.5.pdf).

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
multi-dataset CSV, and combined HTML report described above. Each dataset also has `run_metadata/`
with its normalized sample-sheet row, resolved configuration, command, and
dataset log. Within the harmonisation directory:

- `<dataset>_screen_report.txt` is the readable terminal transcript described
  above;
- `<dataset>_harmonisation_report.html` is the standalone dataset report;
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
The strand section also distinguishes the number of population-reference
unmatched variants detected, retained by explicit opt-in, and removed under the
default policy. Retained rows use the auditable actions
`reference_unmatched_retained` or
`reference_unmatched_retained_reverse_complement`; their population-reference
AF is missing, so they are excluded from study-versus-panel AF comparison. They
remain subject to every later EAF and completeness rule and to the GWAS-to-VCF
genome-FASTA and exact record-count gates. The chromosome table in the HTML
report shows retained and removed unmatched counts separately.

## What to review before downstream use

Confirm final status, inferred input build, complete chromosome set, primary
VCF/index integrity, rejected fractions and reasons, liftover loss, raw and
virtual QC counts, effective-N distribution, study/reference AF concordance,
and any concordance-validation failures. Keep the raw VCF, reports, manifest,
resolved YAML, and logs together.
