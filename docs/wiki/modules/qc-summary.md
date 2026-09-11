# QC Summary

## Purpose

The QC module assesses an existing genotype-free GWAS-VCF and writes raw,
rule-level, and combined virtual-QC metrics. It does not alter the input or
create a filtered VCF.

## What the analysis does

PostGWAS evaluates every configured QC rule independently against every input
GWAS-VCF record, then reports the conjunction of those rules as a virtual
passing subset. It measures P-value, allele-frequency, imputation quality score,
variant-type, palindromic-SNP, MHC, sample-size, and study/reference frequency
evidence without silently discarding or rewriting variants.

Harmonisation calls the same QC-owned assessment internally, so direct and
integrated reports use one scientific policy and one implementation. Reported
rule counts may overlap; the final accounting separately records unique failures
and records matching multiple rules. Each rule also reports its input-wide
percentage, variants excluded only by that rule, and variants that also match
another active rule. Threshold details distinguish values below the configured
minimum from values above the configured maximum.

Rules are classified by scientific purpose. Data completeness or
imputation-quality checks are kept separate from reference concordance, allele
compatibility, and analysis-scope decisions such as SNP-only or MHC-region
exclusion. Consequently, the reports say that matching variants would be
excluded from the virtual subset; they do not claim that every scope exclusion
is a defective variant.

PostGWAS deliberately does not run `bcftools stats` here. That command normally
derives allele-frequency bins and singleton counts from `INFO/AC` and `INFO/AN`,
or from `FORMAT/GT`. Genotype-free GWAS-VCFs normally carry summary-statistic
`FORMAT/AF` instead, so those default bins would not describe the study
frequency. This module reports only metrics derived from explicitly configured
GWAS-VCF fields.

## When to use it

Use QC summary after harmonisation to assess data quality, understand the effect
of candidate filtering rules, or compare raw and virtually passing variants.
Use `postgwas sumstat_filter` when the analysis requires a materialized filtered
VCF. Do not use QC summary to repair an incorrect genome build, allele
orientation, reference population, or sample mapping.

## Input requirements

- A non-empty PostGWAS-harmonised, genotype-free GWAS-VCF whose header declares every configured
  FORMAT and INFO field. An undeclared tag is a structural input error; a `.` in
  a declared field remains an ordinary per-variant missing value and follows the
  configured missing-value action.
- Exactly one supported build declaration in the VCF header (by default,
  `##genome_build=GRCh37` or `##genome_build=GRCh38`). QC infers the build from
  this declaration and uses the matching MHC interval when that rule is active.
- A dataset ID and output directory.
- `bcftools` available through the configured resource command (the packaged
  value is `bcftools`, resolved from `PATH`), including a version response and
  access to the input VCF.
- Exactly one VCF sample column. A lone sample ID that differs from the run
  dataset ID is used with a warning; a zero-sample or multi-sample VCF fails.
- A reference-frequency INFO tag whose allele orientation and population are
  compatible with the study frequency when concordance is assessed.

The public QC command requires unique, non-empty PostGWAS version, dataset-ID,
and status metadata before extraction, resume, or scientific-output replacement.
Their default names are `postgwas_version`, `postgwas_dataset_id`, and
`postgwas_vcf_status`, resolved from the shared
`modules.formatting.input_contract.provenance_headers` configuration. Missing,
blank, duplicated, or malformed origin metadata rejects the study VCF; rerun
PostGWAS harmonisation on the original summary statistics rather than adding
headers manually. Required INFO/FORMAT declarations remain mandatory.
When configured PostGWAS scientific metadata is present, the report captures
its AF source and meaning plus the imputation quality score source,
interpretation, and output field. Missing additional field provenance is
reported explicitly rather than guessed. QC reports only the build of the GWAS-VCF being assessed;
it does not repeat the original summary-statistics build, harmonisation output
build, or liftover history. The harmonisation output-build header is still
checked internally, and QC stops before analysis if it conflicts with the
assessed VCF's `##genome_build` declaration.

## Command

The standalone command is `postgwas qc`; the pipeline target is `qc_summary`.
Inspect the current direct options with:

```console
postgwas qc --help
```

## Minimal example

```console
postgwas qc \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --dataset-id STUDY \
  --output-directory results
```

## Full example

Export and edit the canonical YAML, then supply only the run-specific inputs on
the command line:

```console
postgwas config export \
  --module qc_summary \
  --style full \
  --output qc_summary.yaml

postgwas qc \
  --run-config qc_summary.yaml \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --reference-af-column EUR \
  --minimum-maf 0.01 \
  --minimum-info 0.7 \
  --maximum-info 1.05 \
  --missing-af-action remove \
  --missing-info-action remove \
  --remove-palindromic \
  --remove-mhc \
  --sample-size-reference-quantile 0.9 \
  --sample-size-minimum-fraction 0.6666666667 \
  --dataset-id STUDY \
  --output-directory results
```

Explicit CLI values override the corresponding resolved YAML value. Omitted
CLI options do not replace YAML settings.

The commonly adjusted scientific policy is available directly in both `qc`
and `sumstat_filter`: P-value, MAF, INFO, missing-value actions, reference-AF
concordance, indel handling, palindromic handling, and MHC handling use the same
public option names. QC also exposes its descriptive low-Neff reference
quantile and minimum fraction. `--include-indels`, `--remove-palindromic`, and
`--remove-mhc` are positive presence flags. Their packaged YAML defaults are
SNP-only assessment with no palindromic or MHC removal, so omitting a flag
leaves that behavior unchanged; supplying it opts into the named behavior. To
disable a behavior that was enabled in YAML, edit the run configuration rather
than passing a negative CLI flag. VCF field mappings, MHC coordinates, output paths, executable commands,
and other advanced settings remain YAML-only to keep the direct interface
focused. QC has no `--genome-build` or `--bcftools` option: it infers the build
from the VCF header and resolves bcftools from `resources.executables.bcftools`.
When QC and filtering are selected in the same pipeline, an explicit shared
CLI policy option applies to both modules; use their separate YAML sections
when the two modules should use different values.

## Parameters

`modules.qc_summary` is the canonical schema-validated configuration. Export it
instead of copying defaults from an older run. Its main groups are:

- `inputs`, `reference_af_column`, and `output_directory`;
- `rules`, including P-value, AF, INFO, frequency-difference, indel,
  palindromic, MHC, and effective-sample-size policies;
- `vcf_fields`, which defines the exact bcftools queries;
- `table`, which defines temporary-table parsing; and
- `output_layout`, which defines every report, log, configuration, temporary,
  and completion-manifest path.

Missing-value actions are explicit for P value, AF, and INFO. MHC intervals are
one-based and inclusive in the configuration. The packaged Genome Reference
Consortium regions are GRCh37 `chr6:28,477,797-33,448,354` and GRCh38
`chr6:28,510,120-33,480,577`; change the appropriate build entry only when the
analysis protocol justifies a different interval.

Former harmonisation-only QC/filter keys are not a second source of truth.
Equivalent behavior now belongs under `modules.qc_summary`, including
`rules.maf_min`, `rules.info_min`, `rules.info_max`,
`rules.maximum_af_difference`, the missing-value actions, palindromic policy,
MHC regions, and `rules.sample_size_outlier_standard_deviations`.
The QC assessment and physical filtering module consume one shared scientific
policy contract and are protected by backend-parity tests. Their packaged
upper-INFO defaults remain intentionally different: QC uses `info_max: 1.05`
for a strict report-only virtual subset, while filtering uses `info_max: null`
until the user explicitly requests physical removal. Both schemas accept an
explicit upper limit through 2 for MaCH Rsq data.
Low effective sample size is descriptive and configuration-driven:
`rules.sample_size_reference_quantile` selects the raw usable-Neff reference
quantile using linear interpolation, and
`rules.sample_size_minimum_fraction_of_reference` multiplies that reference to
produce the lower threshold. The packaged values are 0.9 and two-thirds,
respectively. Linear interpolation matches the pandas default used by LDSC.

## Processing steps

1. Resolve and schema-validate the global and QC-summary configuration.
2. Validate the VCF, exact single-sample contract, dataset-ID comparison,
   reference-frequency tag, and the configured bcftools executable.
3. Read the VCF header once, infer exactly one supported genome build from the
   configured build declaration, resolve the build-specific output paths, and
   require every configured INFO/FORMAT query tag to be declared. Report all
   missing declarations together and stop before resume, extraction, or
   scientific-output replacement when the contract fails. A missing,
   duplicate, or unsupported build declaration stops even earlier and is
   recorded in the build-independent preflight log.
   In the same header read, validate the public PostGWAS-origin contract and
   capture the report-relevant PostGWAS AF and score
   provenance using the canonical harmonisation provenance-name mapping.
4. Resolve the configured QC policy into the exact active rule criteria and
   inactive optional rules, and validate that the rule set agrees with the
   shared filtering-policy contract.
5. Validate or restart from the global checksum-backed completion boundary.
6. Run one strict full-record `bcftools query` into an isolated temporary table.
7. Start a fresh aggregation worker with `POLARS_MAX_THREADS` set from the
   resolved `execution.threads` value before Polars is imported, verify the
   worker's actual thread-pool size, and scan the temporary table twice with
   Polars streaming aggregation. The first
   pass calculates the raw and QC-passed effective-sample-size mean and sample
   standard deviation plus the configured raw reference quantile. The second
   uses those scalar values to calculate every independent rule, combined
   pass/fail status, upper-outlier and low-Neff counts, and virtual-subset
   metrics without broadcasting statistics over every row.
8. Reconcile input, passing, uniquely excluded, multiply matched, per-rule
   unique-only, and per-rule overlap counts.
9. Render the metric TSV, rule TSV, assessment JSON, summary CSV, and detailed
   self-contained HTML report from that one reconciled evidence object. The
   report builders do not reopen the VCF or calculate another set of metrics.
10. Publish all five scientific reports as one rollback-safe artifact set, then
   write completion provenance, log records, and the aligned terminal summary.

Like the MAGMA module, a fresh direct QC run shows six durable stages that name
the real work: input GWAS-VCF validation; QC-rule resolution and validation;
conversion of the required VCF fields to a temporary TSV; rule evaluation and
metric calculation; report generation and validation; and rollback-safe report
publication. A failed stage remains below 100% and names the stage that failed.
The progress events and their durations are also retained in the canonical log.
A validated resume uses a separate four-stage resume plan and never claims that
the VCF was reconverted or reassessed.

The terminal and canonical log use the filtering module's presentation
hierarchy. `Input VCF validation` follows stage 1 and identifies the VCF,
header-declared build and contig count, field-contract status, field provenance,
and configured AF and imputation-quality fields. `QC assessment plan ·
conditions evaluated independently` follows stage 2 and uses filtering's plain
numbered action-list style. The final `GWAS-VCF quality-control summary` is shown
only after the complete report set is validated and published. It contains the
same filtering-style sections: scientific results grouped by rule family, exact
accounting, saved reports, and the final outcome. The measured variant total is
reused from the assessment's single VCF extraction, so presentation causes no
additional VCF scan.

Major QC headings use distinct light foreground colours for validation,
policy decisions, scientific results, numerical accounting, saved reports, and
the final outcome. Only headings and semantic field labels are coloured; data
values remain neutral for readability. Canonical log files retain the identical
text hierarchy without terminal colour-control sequences.

The plan and result sections are generated from the resolved policy rather than
a fixed list. Plan entries say `EXCLUDE` because a match is excluded only from
the reported virtual QC-passed subset; the input VCF remains unchanged. Missing-
value decisions are listed separately from numerical conditions, both the lower
and upper imputation-quality limits are explicit, and inactive optional rules
remain outside the numbered plan. Numbers preserve readable policy order only:
unlike filtering's first-failure attribution, QC evaluates every active rule
independently. `QC results by rule` therefore reports each rule's total failures,
input-wide percentage, variants failing only that rule, variants also failing
other rules, and every directional trigger count, including zero. These
overlapping counts must not be added. Low Neff remains a separate outcome
diagnostic because it is not an exclusion rule.

The canonical QC log records the same validation, plan, grouped result,
accounting, report, and final-outcome blocks shown on screen. It then retains a
`Detailed QC audit` containing required tag names, provenance interpretation,
raw distributions, each rule's purpose and criterion, directional triggers,
unique and overlapping impact, virtual-subset metrics, execution details, and
every output. The TSV, CSV, JSON, and detailed HTML reports retain the same
evidence. This keeps monitoring consistent with filtering without discarding
audit information or recalculating a metric.

The collector uses only a streaming parameter explicitly declared by the
installed Polars API: `streaming=True` for the shipped 0.20 release or
`engine="streaming"` for newer releases. It fails instead of passing an unknown
keyword that could be silently ignored. This follows the
[Polars streaming execution contract](https://docs.pola.rs/user-guide/concepts/streaming/).
Polars documents that its pool can be overridden only through
`POLARS_MAX_THREADS` before process start; the requested and observed worker
pool sizes and enforcement result are written to the TSV, JSON, canonical log,
and completion manifest. See the
[Polars thread-pool contract](https://docs.pola.rs/api/python/stable/reference/api/polars.thread_pool_size.html).

The temporary table is removed whether assessment succeeds or fails. For raw
and virtual-QC records, metrics include record and SNP/non-SNP counts,
transition/transversion ratio, AF and INFO missingness, comparable and
discordant frequency pairs, and effective-sample-size availability, range,
mean, sample standard deviation, upper-outlier threshold and count, raw
reference-quantile value, low-Neff threshold, and the count and fraction below
that threshold. The lower threshold is calculated once from raw usable Neff and
reused unchanged for the raw and virtual-QC stages, so their counts are directly
comparable. Low-Neff variants are reported but do not fail an active QC rule or
alter the input VCF. The legacy upper-tail count remains an informational
distribution diagnostic; it does not produce a warning or filter a variant.

## Outputs

With the packaged output layout, the requested output directory contains:

| Output pattern | Description | Important content | Interpretation |
|---|---|---|---|
| `qc_summary/<dataset>_<build>_vcf_qc_metrics.tsv` | Validation, raw, and virtual-subset metrics | Header-contract evidence, counts, missingness, Ti/Tv, AF concordance, effective sample size | Descriptive QC summaries, not filtered data |
| `qc_summary/<dataset>_<build>_vcf_qc_rule_results.tsv` | Per-rule audit | Scientific category and purpose, criterion, decision, total/percentage matched, unique-only impact, overlap, and directional detail counts | Rule totals can overlap |
| `qc_summary/<dataset>_<build>_qc_assessment.json` | Machine-readable assessment | Header contract, PostGWAS scientific provenance, numbered virtual decision plan, active and inactive rules, directional details, complete metrics, and balanced accounting | Canonical reusable result contract |
| `qc_summary/<dataset>_<build>_qc_summary.csv` | Portable combined summary | Overall accounting, input-versus-passing metrics, rule/detail rows, unique/overlap impact, inactive rules, VCF provenance, report paths, and thread provenance | Spreadsheet- and program-friendly view of the same assessment |
| `reports/<dataset>_<build>_qc_report.html` | Detailed self-contained report | Key findings, scientific field provenance, readable per-rule cards, directional triggers, unique/overlap impact, inactive rules, exact accounting, resolved policy and fields, execution provenance, and linked outputs | Human-readable review; no VCF rescan |
| `qc_summary/<dataset>_<build>_qc_resolved.yaml` | Resolved run configuration | Effective global, resource, and module values | Reproduction metadata |
| `qc_summary/<dataset>_qc_preflight.log` | Build-independent failure log | Configuration, executable, input, and VCF-build inference failures that occur before a build-specific log path is safe | Actionable early-failure audit; not written for an ordinary successful run |
| `qc_summary/<dataset>_<build>_qc.log` | Canonical service log | Inputs, parameters, decisions, outputs, and final status | Primary execution audit |
| `qc_summary/<dataset>_<build>_qc_complete.yaml` | Completion manifest | Configuration digest, input/output fingerprints, aggregation strategy, and software versions | Evidence required for validated resume |

No `.stats` file is written beside the input VCF.

## QC and logs

The canonical log records the header-inferred build and its source, resolved rules and fields, VCF metadata,
bcftools path and version, command execution, row accounting, report paths,
warnings, restart decisions, and final completion status. Header validation
records every required tag, declared INFO/FORMAT counts, and either `PASSED` or
all missing declarations. A failed contract also records the run as `FAILED`;
it never records `COMPLETED` or publishes new assessment reports. Successful
terminal and tabular summaries group the input filename, build, declared
contigs, measured variant total, header-contract result, required declared
fields, and available AF/score provenance into one validated-input section. The
log separately records every active rule's category, criterion,
total match fraction, unique-only count and overlap count, and every inactive
rule with its reason. Screen fields share one aligned value column and
the final section names the CSV and HTML reports. The completion manifest
fingerprints the input VCF and all five scientific assessment reports, and
records the PostGWAS, Python,
Polars, and bcftools versions plus the streaming aggregation strategy.
The result log also records the resolved low-Neff quantile and fraction, the
resulting reference and threshold, and raw and virtual-QC low-Neff counts and
fractions. The completion manifest records the threshold and both counts.

Global resume is enabled by default. A checksum-valid completed result is reused.
A real partial result or changed parameter/input warns and restarts from the
safe QC-summary boundary, replacing only checksum-matching PostGWAS-owned files.
An externally modified output is preserved and refused instead of being deleted
automatically. `--overwrite` remains the explicit forced replacement.

## Interpretation

Reference-frequency concordance is meaningful only when study and reference
frequencies use the same effect/alternate-allele orientation, genome build, and
compatible population. “QC-passed” describes the combined virtual assessment,
not a new VCF and not proof that the study is scientifically suitable for every
downstream method.

The imputation quality score field is not automatically a study-measured
statistic. PostGWAS harmonisation may populate it from an external reference
proxy or a user-assigned fixed value. The generated report therefore uses
“imputation quality score” while displaying the VCF's `info_source`,
`info_interpretation`, and output-field provenance. Standard INFO and MaCH Rsq
interpretation must follow that recorded provenance and the resolved thresholds;
the QC report does not infer score type from its numerical range.

A low-Neff warning means that a variant was analysed in a substantially smaller
effective sample than the dataset's well-covered variants. This can reflect
legitimate variant availability, cohort participation, chromosome-specific
sample size, or genotyping/imputation loss; it is not automatically a data
error. Review the affected fraction and downstream method before configuring a
separate filtering module to remove those variants. The packaged threshold
follows the low-sample-size convention used by LDSC: two-thirds of the raw 90th
percentile. PostGWAS reports it rather than silently applying LDSC's removal
policy to every downstream analysis.

## Common problems

| Problem or error | Likely cause | How to check | Solution |
|---|---|---|---|
| Required VCF declaration is absent | The input schema, configured `vcf_fields`, or `reference_af_column` does not match | Find `VALIDATION vcf_qc_header_contract` in the QC log; every missing INFO/FORMAT tag is reported together | Correct the VCF header/input or explicitly configure the real field; do not treat an undeclared tag as ordinary `.` missingness |
| Genome build cannot be inferred | The configured build declaration is absent, duplicated, or names an unsupported build | Read `qc_summary/<dataset>_qc_preflight.log` and inspect the VCF header | Supply a PostGWAS harmonisation VCF containing exactly one supported `##genome_build=<build>` declaration; do not relabel coordinates at QC time |
| VCF sample contract fails | The VCF has zero or multiple sample columns | Run `bcftools query -l` | Supply a VCF containing exactly one GWAS sample; a lone sample-ID/dataset-ID mismatch is reported as a warning |
| Frequency discordance is high | Allele orientation, build, ancestry, or reference release differs | Compare VCF metadata and reference provenance | Use a scientifically compatible reference; do not relax the threshold blindly |
| Resume refuses an output | A recorded output changed after validation | Read `checkpoint_events.log` and compare the recorded checksum | Preserve/review the file; use `--overwrite` only for intentional replacement |

## Limitations

QC summary assesses records but does not materialize a passing VCF, repair
alleles or coordinates, infer reference compatibility, or establish causality.
Independent rule counts overlap by design, although the reports expose each
rule's unique-only and overlapping components. The public origin gate checks
provenance declarations, not cryptographic authentication or complete scientific
QC. The reusable internal assessment called by harmonisation remains separate
from this public-input gate, and external reference-panel VCFs are not required
to carry study-origin metadata. The VCF is queried only once, while
the extracted temporary table is scanned twice to avoid materialising repeated
full-length statistic columns. The table still requires disk space proportional
to the queried records. The QC log and reports record the Polars version,
streaming collection API, aggregation strategy, and temporary-table scan count.

## Scientific references

- [bcftools query and statistics documentation](https://samtools.github.io/bcftools/bcftools.html)
- [Genome Reference Consortium GRCh37 MHC region](https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37)
- [Genome Reference Consortium GRCh38 MHC region](https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38)
- [GWAS-VCF specification](https://pmc.ncbi.nlm.nih.gov/articles/PMC7805039/)
- [LDSC `munge_sumstats.py` low-sample-size filter](https://github.com/bulik/ldsc/blob/master/munge_sumstats.py)
- [Winkler et al. GWAS meta-analysis quality-control protocol](https://doi.org/10.1038/nprot.2014.071)
