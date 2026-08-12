# Harmonisation configuration

## Canonical command

```console
postgwas harmonisation \
  --sample-sheet studies.csv \
  --run-config postgwas.yaml
```

Explicit CLI values override the user YAML, which overrides packaged YAML.
Argparse does not own configurable defaults.

## Sample sheet

The version-2 sample sheet contains dataset-specific paths, identifiers, and
source-column mappings. `resource_directory` and `output_directory` are run
settings and do not belong in dataset rows.

`external_eaf_file` and `external_info_file` accept either one existing table
containing all chromosomes or an explicit path template containing
`{chromosome}` and optionally `{build}`. A common filename prefix alone is not
resolved implicitly, because selecting a reference file by guesswork could use
the wrong chromosome or panel.

`trait_type`, `effect_type`, `p_value_type`, and `delimiter` are inferable.
Missing values become `auto`. Invalid values produce a normalization warning
and become `auto`. Scientific inference and comparison with explicit values is
a later preflight stage.

Explicit non-`auto` sample-sheet declarations have final precedence over the
generic YAML values for `effect.type`, `pvalue.type`, and `input.delimiter`.
PostGWAS still runs automatic inspection as an independent cross-check. When
the two disagree, `validation.declaration_mismatch_action` controls the result:
the default `warn` prints a bright-red terminal warning, records the mismatch,
and continues with the declaration; `fail` stops the dataset before chromosome
processing. Automatic inspection is evidence, not permission to silently
replace an explicit scientific declaration.

### Strand consensus and per-variant orientation

Dataset step 05 tests each candidate genome build with forward, swapped,
reverse-complement and reverse-complement-swapped allele matches. The same
coordinate joins provide the non-palindromic forward/reverse counts used at
dataset step 06, so the build references are not read a second time. A
consensus requires `strand.min_informative_variants` and the dominant fraction
configured by `strand.consensus_threshold`. Palindromic SNPs never contribute
because their allele letters cannot distinguish the competing strands.

Automatic build inference requires all three independent evidence checks: at
least `build.min_match_count` allele-compatible rows (default 100), at least
`build.min_match_fraction` of coordinate-testable rows matching the winning
build (default 0.80), and at least `build.confidence_ratio` of all allele
matches supporting that build (default 0.90). A row is coordinate-testable only
when its exact chromosome and position occurs in at least one configured build
reference; rows absent from both references are reported separately and do not
dilute the fraction. Each study row counts once per build despite duplicate or
multiallelic reference records, and exact duplicate reference rows are removed
by default before joining. An explicit `build.mode` remains authoritative and
bypasses the automatic evidence thresholds while still collecting strand
evidence from the declared build reference.

At chromosome step 04, the supplied raw frequency reference remains an
ordinary `CHROM`, `POS`, `REF`, `ALT`, population-AF table. PostGWAS joins it
once to the study on chromosome and position, then uses vectorised Polars
expressions to classify direct and swapped matches for every variant and the
two reverse-complement matches for SNVs. Reverse-complement matching is not
applied to indels because it is not a safe substitute for VCF normalization.
It does not create or require a four-times-expanded reference file. The
selected result is recorded in
`strand_action`, with aligned reference ALT frequency in
`strand_reference_af`.

Chromosome step 03 has already converted an odds ratio and any declared raw
OR-scale SE to canonical log-odds BETA and SE. Complementing both alleles then
preserves BETA, Z and true EAF. Changing the effect allele negates BETA and Z
and replaces true EAF with `1-EAF`; log-scale SE, p-value, INFO and sample size
remain unchanged. A confirmed MAF is never inverted because it is not tied to
the listed effect allele. If the second-level reference comparison instead
proves that a MAF-like column is true EAF, swapped rows are inverted at that
point. After chromosome processing, PostGWAS consolidates every completed
chromosome's `maf_reference_decision`. A confirmed MAF or an inconclusive
comparison stops that chromosome; PostGWAS does not guess, export an unoriented
frequency, or apply a deferred frequency flip without evidence. Only unanimous
conclusive EAF evidence produces a successful final dataset frequency type.
The manifest preserves the initial MAF suspicion separately and records the
final type, its source, and the chromosome decision counts. A palindromic A/T
or C/G SNP is oriented only when the full-study non-palindromic evidence has a
strong forward or reverse consensus. Mixed, unresolved, or insufficient
consensus rejects it as `palindromic_ambiguous`; population AF never chooses
its strand. The defaults reject unmatched or ambiguous rows with explicit
provenance. After the final internal or external EAF is known,
`strand_af_difference` compares it with the aligned population AF. A large
difference is supporting QC evidence, not a reason to reorient the row.
Non-palindromic discordance warns by default through
`strand.af_discordance_action`; palindromic discordance rejects by default
through `strand.palindromic_af_discordance_action`. These decisions follow the
[GWAS Catalog summary-statistics harmonisation method](https://www.ebi.ac.uk/gwas/docs/methods/summary-statistics).

External EAF and INFO annotations use the same row-preserving direct/swapped
allele matcher after strand orientation. It normalizes the four join keys,
gives direct matches precedence, performs one Polars left join and preserves
the study row count and order. On a swapped match, EAF becomes
`1-AF`; INFO is allele-independent and remains unchanged. In the normal
sample-sheet workflow the EAF source is exclusive, so an external EAF does not
perform a second frequency-driven strand decision after `strand_action` has
been fixed. INFO retains its independent range and missing-value policies. The
strand, external EAF, and external INFO references use the single
`external_reference` duplicate policy. The defaults keep one row only when the
normalized chromosome, position, alleles and selected numeric annotation are
identical (`exact_duplicate_action: keep_one`). If values differ—including a
missing value beside a finite value—all reference rows for that allele key are
discarded (`non_identical_duplicate_action: discard_all`), so row order never
chooses scientific data. The affected study row is then unmatched and follows
the existing EAF or INFO missing-value policy. Duplicate-group and removed-row
counts are recorded in the chromosome log.

### Effect scale before Z-based recovery

The study-wide effect decision is passed unchanged to every chromosome. At
chromosome step 03, an odds ratio is normalized through the existing
`harmonise_effect_estimates()` implementation, producing the canonical
log-odds beta and applying the configured non-positive-OR and SE-scale rules.
Only then does step 06 derive a missing standard error from `SE = beta / Z`.
Thus an OR uses `SE = ln(OR) / Z`, never `OR / Z`, and the canonical column
mapping prevents a second logarithm later in the workflow. If the input has no
effect column, step 05 is skipped and step 06 creates beta directly from Z.
This follows the standard logistic-regression convention that the reported
coefficient is the log odds ratio and its Z statistic is coefficient divided
by standard error ([Stata logistic-regression interpretation](https://www.stata.com/links/stata-basics/logistic-regression-1-introduction/)).

When `effect.type: auto`, PostGWAS uses only finite numeric effects. An odds-ratio
decision requires both a median in the configured near-one interval (default
`[0.80, 1.25]`) and at most 0.5% combined zero or negative values. A low
non-positive fraction with a median outside that interval is ambiguous; a high
fraction with a median inside it is conflicting. Both conditions stop before
chromosome processing and request an explicit effect type. A high fraction and
a median outside the interval is classified as beta. For a confirmed OR
column, individual non-positive rows are rejected before `ln(OR)` by default.
Whenever this automatic detector supplies the applied study-level answer,
PostGWAS prints and records a prominent warning with the inferred type, median,
non-positive fraction and expected transformation. The warning explains that
values alone cannot distinguish every all-positive beta distribution from an
odds-ratio distribution and recommends an explicit sample-sheet `effect_type`.

### Effect-statistic concordance

Chromosome step 09 handles missing or non-finite BETA, unsafe SE values and
missing or non-finite Z before calculating `Z = BETA / SE`. In the normal
pipeline, step 10 verifies that it received the exact resulting DataFrame under
the same policies and reuses those basic results instead of scanning the three
columns again. Any changed DataFrame, column mapping or policy makes step 10 run
the full checks. The positive, configurable SE division floor used by step 09
is stricter than the scientific requirement `SE > 0`, so reuse cannot admit a
zero, negative or numerically unsafe denominator.

Chromosome step 10 validates the standard Wald relationship `Z = BETA / SE`
when all three statistics are present. The default combined absolute and
relative tolerances allow ordinary publication rounding; discordance is warned
and recorded rather than silently removed. The same step now also enables the
two-sided Z-versus-p-value comparison by default, allowing a factor-of-ten
difference under `validation.z_pval_tolerance_log10: 1.0`. Either check can be
configured as `off`, `warn`, `reject`, or `fail` in the canonical YAML. These
checks follow the summary-statistic consistency approach documented by
[MungeSumstats](https://www.bioconductor.org/packages/release/bioc/vignettes/MungeSumstats/inst/doc/MungeSumstats.html).

Every dataset must provide exactly one study EAF source: an internal column XOR
an external file-plus-column pair. INFO is resolved independently and in a
fixed order: an internal INFO column has first priority, an external INFO file
plus column has second priority, and the explicit command-line fallback
`--fixed-info VALUE` is used only when neither is present. If all three are
absent, preflight fails before chromosome processing. An internal INFO column
also wins when the same sample-sheet row lists an external INFO source; the
ignored source is recorded.

The command-line parser requires `0 <= --fixed-info <= 1`. A constant such as
`--fixed-info 0.99` is assigned to every retained variant in only the affected
dataset and is exported through the ordinary INFO/SI mapping. The dataset log,
chromosome QC and run manifest identify it as a user-assigned constant rather
than a measured per-variant imputation-quality score. It is never activated by
omission and cannot be set through the run YAML. This explicit provenance is
important because INFO is encouraged rather than mandatory in GWAS-SSF and SI
is optional in GWAS-VCF; a constant therefore supplies a requested output value
but does not reconstruct the unavailable variant-level measurement
([GWAS Catalog format](https://www.ebi.ac.uk/gwas/docs/summary-statistics-format),
[GWAS-VCF specification](https://github.com/MRCIEU/gwas-vcf-specification)).

### Sample-size semantics

The sample-sheet `trait_type` is passed to sample-size harmonisation:
`case_control` selects the `binary` policy, `quantitative` selects the
quantitative policy, and `auto` retains the canonical YAML policy. An explicit
non-`auto` dataset declaration has final precedence over a generic policy.

Case-control studies require both case and control counts and use
`4/(1/Ncase + 1/Ncontrol)`. A case count alone fails because it
cannot determine case-control effective sample size. Quantitative studies keep
the existing interface: `control_count_column` or `control_count` represents
the total analysed N and is copied to Neff; supplying any case-count input for
an explicitly quantitative study fails in both sample-sheet validation and the
internal sample-size step. Every configured sample size and Neff must be at
least one. The effective-N convention follows
[Mallard et al.](https://onlinelibrary.wiley.com/doi/full/10.1002/gepi.22609).

At the current sample-sheet preflight stage, referenced data files are never
opened or decompressed. Validation uses filesystem metadata only: the path must
exist, identify a regular file, and have nonzero size.

## Execution defaults

Packaged automatic rules reserve one logical CPU and ten percent of physical
memory. The actual resolved values are shown in `--help` and written to
`resolved_config.yaml`.

The shared `--threads` and `--memory-gb` options are identical for standalone
modules and pipeline mode. Their implementation lives in one common compute
CLI component. No deprecated harmonisation interface is available.

## Runtime prerequisites

Harmonisation requires `bcftools`, `tabix`, and the bcftools `liftover` plugin.
The project container builds matching HTSlib/bcftools 1.23.1 sources and pins the
plugin source to an immutable revision. A local installation can be checked with:

```console
bcftools --version
bcftools plugin -l | grep '^liftover$'
```

After creating and activating the conda environment, install the pinned plugin
with `tools/setup/install_bcftools_liftover.sh`. Optional first and second
arguments override the bcftools version and plugin source revision.

The command verifies these requirements before reading a summary-statistics
file. A missing program or plugin therefore produces an actionable preflight
error in `run_metadata/postgwas.log`, rather than failing after chromosomes
have already been processed.

## Logs and run metadata

Every run writes `run_metadata/postgwas.log` and `resolved_config.yaml` at the
output root. It also writes the configured `runtime.run_summary_file` (packaged
default `harmonisation_run_summary.csv`) with one row per selected dataset.
Its columns follow pipeline order. This run-level CSV records final status;
input-QC-ready SNP and indel/other counts; resolved and detected build/effect/P
types, with explicit automatic-detection flags; the harmonised effect scale;
combined completed-chromosome strand counts and reference files; final-VCF SNP
and indel/other totals; SNP QC pass/fail counts and percentages;
effect-allele-frequency concordance evidence; selected statistical-quality
counts; and the dataset manifest path. It is initialized before dataset
analysis and updated atomically after every dataset, including failure and
interruption states. The combined values reuse manifest counters and the
persisted QC assessment JSON, and perform no additional GWAS or VCF scan. A
command-level preflight failure after the sample sheet has been parsed records
every selected dataset as `PREFLIGHT_FAILED` with the shared reason; property
and chromosome fields remain blank because no scientific analysis occurred.

Normal terminal output is always copied incrementally to the configured
`output_layout.screen_report`. With the packaged layout, each dataset receives
`<output>/<dataset_id>/harmonisation/<dataset_id>_screen_report.txt`. Shared
run preparation and final summary blocks are copied to every selected dataset;
dataset processing blocks are written only to the dataset they describe.
`runtime.display_screen` defaults to `true`. Use `--hide-screen` to suppress
normal progress on stdout while continuing to write the reports, or
`--show-screen` to override a run configuration that disabled display. Fatal
command errors remain visible on stderr.

Each dataset also receives its own run metadata and complete dataset/chromosome
logs. Success, permanent failure, configuration failure, and user interruption
all append an explicit final record; the screen remains limited to progress,
warnings, failures, and result locations.

Harmonisation chromosome logs use the same compact audit vocabulary throughout:

- `STEP`, `INPUT`, and `PARAM` identify the operation, incoming row count, and
  resolved values actually consulted by that step.
- `OBSERVED`, `DECIDE`, `ACTION`, and `RESULT` show the evidence, selected
  analysis path, operation performed, and its measured outcome.
- `PASS` records a completed check that affected no values in one line.
- `OUTPUT` identifies generated artifacts, while `STATUS`, `WARNING`, and
  `FAILED` give the final state and actionable problems.

Long parameter explanations are available from the canonical configuration
export and are not repeated for every chromosome. External-tool transcripts
are kept apart from scientific step logs. GWAS-to-VCF transcripts are written
under each dataset's `harmonisation/logs/adapters/gwas2vcf/` directory;
each contains the exact command, complete tool output, and exit code.

## Internal module boundary

Harmonisation analysis files contain the scientific operation and its public
entry point. Repeated mechanics—policy resolution, null logger/step contexts,
row rejection, optional-value handling, Polars normalization, and checked
commands—are shared rather than reimplemented in every chromosome step.
Harmonisation-specific mechanics live under `harmonisation/shared/`; utilities
that are safe for every PostGWAS module live under `postgwas/core/`. The
GWAS-to-VCF adapter remains isolated under `harmonisation/adapters/`.

## Duplicate variants

Duplicate handling is performed once during dataset-level input validation,
before chromosome partitioning. `duplicates.key` identifies groups; its default
is chromosome, position, effect allele, and other allele. Allele letter case is
ignored, but effect/other allele order is retained because effect type and
reference orientation have not yet been resolved at this stage.

For every duplicate group, PostGWAS compares the non-empty columns selected by
`duplicates.consistency_fields`. If those scientific values agree, one row is
selected using `duplicates.selection_order`: completeness across configured
quality fields, larger per-variant effective sample size, higher INFO, and
finally the earlier physical input row. Case-control effective sample size uses
`4/(1/Ncase + 1/Ncontrol)`; a quantitative-trait N is used directly. The
formula is the usual balanced case-control equivalent sample size
([Mallard et al.](https://onlinelibrary.wiley.com/doi/full/10.1002/gepi.22609)).
The input-order tie-breaker makes the result reproducible without selecting on
statistical significance.

If a group contains different non-empty association values, every row in that
group is rejected as `conflicting_duplicate`; `duplicates.conflicting_action`
can instead stop the dataset after the reports are written. The duplicate TSV
contains every member of each group plus its classification, conflicting
fields, and `duplicate_action`. The row retained for analysis is labelled
`kept`; other consistent rows are labelled `removed`, while conflicting rows
carry the configured conflict action. The rejected-variant file separately
records every row that was actually removed. `duplicate_input_row` is the
one-based physical source row used by concordance to identify that exact row;
concordance does not repeat the duplicate ranking.

Rejected-variant provenance uses the same stable one-based source-row ID. The
input reader captures an immutable snapshot immediately after parsing and before
coordinate, allele, frequency, effect, or standard-error transformations. Only
the snapshot for the active chromosome is loaded when rejected rows are written,
so the rejected file contains the original parsed study values even when the
working row was subsequently swapped, inverted, rescaled, or converted. The
internal chromosome snapshots are compressed Parquet intermediates and are
removed during successful final cleanup.

This policy deliberately does not keep the smallest p-value from conflicting
records, because selecting on significance can bias the retained result. It
supports the one-record-per-variant design of GWAS-VCF for duplicate groups
that share the configured ordered allele key, while retaining an auditable
record of exclusions. Supporting conventions and comparison
implementations are described by the
[GWAS-VCF publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC7805039/), the
[GWAS Catalog summary-statistics format](https://www.ebi.ac.uk/gwas/docs/summary-statistics-format),
and the
[MungeSumstats duplicate checks](https://www.bioconductor.org/packages/release/bioc/vignettes/MungeSumstats/inst/doc/MungeSumstats.html).

## Merged-VCF QC assessment

Before QC, each required merged VCF group is reconciled against the chromosomes
that completed harmonisation. The input-build annotated VCF, target-build lifted
VCF, and raw gwas2vcf VCF are declared by
`vcf_processing.required_merge_groups` and are all required. A missing or invalid chromosome
input, failed concatenation, missing or undersized merged output, missing tabix
index, or unreadable indexed record count is recorded as a merge failure.
`vcf.on_merge_failure` defaults to `fail`, so incomplete genome-wide outputs do
not continue into QC. Per-chromosome inputs are retained for retry. The merged
not-lifted VCF remains optional when no chromosome produced rejected variants.

The Python GWAS-to-VCF adapter receives the study table and source-build FASTA,
but not dbSNP. It preserves the configured study variant identifier independently
in `FORMAT/ID`; arbitrary non-empty identifiers are accepted and made VCF-safe.
After source-build normalization, bcftools removes the temporary record ID and
assigns the record-level dbSNP ID using exact normalized chromosome, position,
REF and ALT matching. This prevents position-only assignment at multiallelic
sites. IDs still missing after the dbSNP match are filled with the configured
`vcf_processing.missing_id_format`. Its default,
`+%CHROM\_%POS\_%REF\_%ALT`, follows the bcftools `--set-id` syntax and escapes
literal underscores so they are not parsed as part of field names. The matching
and replacement behavior follows the official
[bcftools annotate documentation](https://samtools.github.io/bcftools/bcftools#annotate).

Harmonisation writes the raw merged VCF and does not create a second,
QC-filtered VCF. Immediately after input-build concatenation, the optional
population-frequency check runs a narrow `bcftools query` against that
unfiltered VCF. It extracts only the five configured INFO frequencies: the
study `INFO/AF` plus AFR, EAS, EUR and SAS reference AF. No coordinate or allele
column is extracted. Missingness is reported for every field. Correlations and
mean absolute differences use only the common records where all five
frequencies are finite and between zero and one, so different population-tag
missingness cannot give each population a different comparison set.

The closest label is reported only when the configured minimum comparable
count, minimum correlation, and best-versus-second correlation gap pass, and
the highest-correlation population agrees with the smallest absolute AF
difference when `require_mae_agreement` is enabled. This is descriptive
reference-frequency compatibility QC, not individual-level or cohort-level
ancestry inference; it never changes a population label or filters a variant.
The four packaged population labels are a configured subset of the 1000
Genomes Phase 3 super-populations described by the
[1000 Genomes Project](https://doi.org/10.1038/nature15393).

If the sample sheet supplies `external_eaf_file` or `external_info_file`, only
its basename is normalized to uppercase alphanumeric tokens and inspected for
an exact configured population token. A strong EUR similarity result paired
with an external filename that contains the exact token `AFR`, for example,
produces a warning and continues; `AFRICA` does not count as `AFR`. File
contents, INFO-score values and path directories are not used for this filename
check. No recognizable token means no filename comparison. The terminal and
JSON report state whether the selected comparison-AF column and each recognized
external filename match the closest population. The structured result is
written to `<dataset>_<build>_population_frequency_qc.json`.

After output validation, a six-card final QC takeaway reuses the already
calculated dataset, chromosome, merged-VCF and virtual-QC results; it does not
query or scan the VCF again. Variant accounting reports input rows, rows read,
rows sent to GWAS-to-VCF, final merged-VCF records, SNP/indel content, and
QC-passed and failed SNPs out of all raw SNPs. Reference-AF concordance uses
only variants with both the study AF and selected population AF present. It
reports concordant and mismatched counts and percentages to four decimal
places, the configured absolute-AF
difference cutoff, and study/reference missingness separately. The remaining
takeaways summarize build and strand evidence, effect/Neff/INFO quality,
rejection provenance, explicit count reconciliation, and separate chromosome,
merged-VCF and output integrity. “QC-passed” is a virtual assessment; the raw
merged VCF remains unchanged.

The existing merged-VCF rule assessment then runs its own `bcftools query` data
pass against the build selected by `qc.target_build`, extracting only the fields
required for that assessment into a temporary TSV. A streaming Polars
aggregation performs four related calculations without rereading its VCF:

1. summary-statistic flow before VCF creation;
2. raw merged-VCF metrics, including variant type, Ts/Tv, missing QC fields and
   study/reference allele-frequency concordance;
3. every active configured QC rule evaluated independently against all raw VCF
   records, so rule-specific counts remain scientifically interpretable; and
4. metrics for the final virtual QC-passed subset, created by applying all
   active rule masks together.

The same extraction includes the configured effective-sample-size FORMAT field
(default `FORMAT/NEF`). For both the raw records and the final virtual subset,
PostGWAS reports usable and missing/invalid values, minimum, maximum, mean,
sample standard deviation, and the number above that stage's mean plus
`qc.sample_size_outlier_standard_deviations` standard deviations. The default is
5, matching
[MungeSumstats `N_std`](https://www.bioconductor.org/packages/release/bioc/manuals/MungeSumstats/man/MungeSumstats.pdf).
This is descriptive QC: unusually large Neff values are reported but are not
silently removed.

“QC-passed” means that a raw-VCF record passed every active rule under
`policies.filter`. It does not identify another output VCF. The terminal report
shows how many raw records fail each independent condition and explicitly notes
that these counts can overlap. Only the final combined mask determines the
excluded and QC-passed totals.

Persistent reports are written to the dataset's `qc_summary/` directory:

- `<dataset>_chromosomewise_harmonisation_metrics.tsv` contains flattened
  chromosome-wise and dataset-level harmonisation metrics;
- `<dataset>_pre_vcf_column_statistics.tsv` contains per-chromosome statistics
  for the harmonised columns exported for VCF creation;
- `<dataset>_gwas2vcf_column_mapping.json` contains the validated common
  mapping from GWAS-to-VCF field names to exported column positions;
- `<dataset>_<build>_vcf_qc_metrics.tsv` contains raw and final metrics;
- `<dataset>_<build>_vcf_qc_rule_results.tsv` contains raw-VCF rule criteria,
  actions, and affected-variant counts;
- `<dataset>_<build>_qc_assessment.json` contains the complete structured result.

The temporary extracted TSV is always removed, including when assessment fails.
The raw merged VCF is unchanged and remains the VCF returned by harmonisation.
The bcftools query expressions are declared under `vcf_processing.qc_fields`,
and the three persistent report paths plus the temporary-file prefix are
declared under `output_layout`. Screen labels are derived from the resolved VCF
mapping, so a configured tag change is reflected in both analysis and reporting.

The supported comparison-frequency panel names and their resource examples are
declared under `comparison_af.available_sources` and
`comparison_af.resource_examples`; the CLI choices and help are generated from
those values rather than maintained independently. This indexed VCF supplies
population INFO annotations and post-merge frequency QC. The packaged default
is ALFA with the EUR population field.

When the study-wide frequency distribution is MAF-like, chromosome step 04
checks it against `default_eaf.column` from the tabular panel selected by
`default_eaf.source` and declared by `resource_layout.default_eaf`. The shipped
default is ALFA EUR. The validated tabular choices are `ALFA`, `wgs_ukb`,
`panukb`, `1000G`, and `fingen`; their exact resource examples are maintained
under `default_eaf.available_sources` and `default_eaf.resource_examples` in
the canonical YAML. Every choice uses the same configured chromosome TSV path
contract and must contain the selected population column. This role remains
independent of the comparison VCF annotation panel. Structural columns and the
tab delimiter come from `default_eaf_mapping`. Autosomes are available for all
five configured panels and both supported builds; sex-chromosome coverage is
panel- and build-dependent, so dataset preflight verifies every chromosome
actually required by the study. Direct allele matches use the ALT frequency;
swapped matches use `1-AF`. The study column is
accepted as EAF only after at least `eaf.maf_reference_min_overlap`
non-palindromic matches, when no more than
`eaf.reference_minor_fraction_cutoff` of the aligned reference EAF values are
at or below 0.5, and when its mean absolute error as EAF is lower than its
error as folded MAF by more than
`eaf.maf_reference_error_margin`. The reverse error separation confirms MAF.
Too few matches, a reference set dominated by minor effect alleles, or errors
closer than the required margin are inconclusive. Confirmed MAF and
inconclusive results both stop with the measured evidence; neither is silently
used as effect-aligned EAF. Palindromic A/T and C/G variants do not contribute
to this decision.

All non-VCF reference tables use the same shared, compressed-file-aware reader.
`external_eaf_mapping`, `external_info_mapping`, and `build_check_mapping`
declare their chromosome, position, first allele, second allele, and delimiter
in the canonical YAML. A delimiter can be explicit or `auto`; automatic
detection is limited to the configured `input.delimiter_candidates` and fails
when no candidate produces a consistent table. The genome-build check paths are
expanded from `resource_layout.build_check` for every build named by
`vcf_processing.target_builds`. Build inference therefore uses configured build
names and reference schemas, and counts as comparable only study chromosomes
that occur in at least one configured build reference. It accepts direct and
swapped allele order but counts each study row at most once for each candidate
build; repeated reference rows or both reference orientations cannot inflate
the evidence. This follows the GWAS Catalog convention of one row per variant
while allowing either other/effect-allele order ([format specification](https://www.ebi.ac.uk/gwas/docs/summary-statistics-format),
[harmonisation methods](https://www.ebi.ac.uk/gwas/docs/methods/summary-statistics)).

Once build inference and study-wide decisions are complete, harmonisation
expands `resource_layout` for only the chromosomes present in that dataset. A
mandatory preflight then validates all resolved files together: VCF/FASTA
indexes and contig compatibility, configured population INFO definitions,
reference-table schemas and data presence, GFF structure, and source-to-target
chain mappings. Shared resources are read once. Any problem fails the dataset
before chromosome partitions or workers are created, and the complete grouped
problem list is recorded in the dataset log. No separate resource-preflight
policy is exposed because these checks protect inputs required unconditionally
by the active bcftools/GWAS-to-VCF commands.

## Input-to-VCF accuracy validation

Validation is off by default. Add `--validate` to a harmonisation run to audit
every dataset after its same-build, unfiltered merged VCF has been created:

```console
postgwas harmonisation \
  --sample-sheet studies.csv \
  --run-config harmonisation.yaml \
  --output-directory results \
  --validate
```

Use `--dataset-id` to limit that command to one sample-sheet row. An existing
harmonisation result can be audited independently:

```console
postgwas --validate \
  --sample-sheet studies.csv \
  --dataset-id study1 \
  --vcf study1_GRCh37_merged.vcf.gz \
  --run-config harmonisation.yaml \
  --output-directory results
```

The audit matches chromosome, position, and alleles in direct, swapped, strand
complement, and complement-plus-swapped orientations. The allele-aware match
total includes all of those recognized orientations, while each orientation is
also counted explicitly. SNPs, indels, and equal-length multi-base variants are
reported separately. Input-only and VCF-only variants are descriptive evidence:
their number never fails this value-preservation audit because VCF normalization,
multiallelic splitting, and earlier filtering can change record representation.

The screen and TSV summary first report input composition, final-VCF
composition, and the allele-aware comparison. A single common-variant
percentage is calculated as the intersection divided by the unique
chromosome-position-allele union. Input-only and VCF-only percentages use their
own source totals. Match orientation and matched-value concordance remain
separate subsections because a correctly identified swapped allele pair can
still contain a wrongly signed effect or uninverted frequency.

After allele-aware matching, only unmatched records enter position diagnostics.
They are joined by chromosome and position. Only positions containing exactly
one unmatched input record and one unmatched VCF record of the same variant type
are compared; one-to-one type mismatches and multiallelic positions are labelled
ambiguous rather than paired arbitrarily. Because allele
orientation is unknown at this stage, position-level effect and Z comparisons
use absolute magnitudes, allele frequency is folded to minor frequency, and SE
is compared directly. These checks are reported and never establish variant
identity or fail the dataset. Truly input-specific and VCF-specific variant and
position counts remain visible. This restrained diagnostic is useful because
[`bcftools norm -f`](https://samtools.github.io/bcftools/bcftools#norm)
left-aligns and normalizes indels before concordance.

For allele-aware matches, PostGWAS compares log-scale effect, canonical
log-scale standard error, effect-allele frequency, and Z score after orienting
allele-specific values to the VCF ALT allele. Odds ratios are converted to
log(OR); a study-level `as_given` OR-scale SE is divided by OR before comparison.
Swapped matches negate BETA and Z and use `1-EAF`; log-scale SE remains
unchanged. When the input has no Z column, Z is calculated from effect/standard
error, or from signed effect and the configured p-value scale when standard
error is unavailable. When input SE is absent but BETA and a non-zero Z are
available, the expected canonical SE is reconstructed as `abs(BETA/Z)` for the
comparison. Global metrics and separate SNP, indel, and other-variant metrics
are written. The default `palindromic_action: compare_resolved` includes every
retained A/T or C/G SNP using the same strong forward/reverse study consensus
recorded in the run manifest; it reconstructs aligned versus swapped direction
without using population AF as a strand decision. If that provenance is absent,
the validator does not guess and excludes only orientation-sensitive values.
`exclude` always skips those values, while `compare_as_listed` is an explicit
override for independently pre-aligned input and does not use harmonisation's
strand proof.
Palindromic comparisons and all numerical tolerances are controlled by
`concordance_validation` in the harmonisation YAML. Input parsing and p-value
interpretation reuse the canonical `input` and `pvalue` policies rather than
defining a second set of values.
Concordance reads the final reference-resolved frequency type from the run
manifest. It performs effect-allele-oriented EAF comparison after reference
confirmation of EAF, preventing `0.2` and an incorrectly inverted `0.8` from
appearing concordant after MAF folding. If the final type is unresolved, only
the allele-frequency comparison is skipped with a warning; variant matching,
effect, and Z-score checks still run.

Reports are written below
`<output>/<dataset>/harmonisation/qc_summary/concordance/`.
The status summary, matched-value mismatches, input-only variants, VCF-only
variants, duplicate VCF records, and unambiguous same-position diagnostic pairs
are always written; the full allele-aware matched table is optional. Unmatched
records produce `WARNING`, not `FAIL`. A failure is reserved for configured
matched-value disagreement, duplicate VCF keys, invalid VCF records, or an
execution/preflight error. A complete audit log is written to the sibling
`logs/` directory even when preflight or analysis fails.

Terminal reporting follows the scientific order: (1) input composition, (2)
final-VCF composition, (3) allele-aware common/orientation/value/unmatched
results, (4) position diagnostics restricted to allele-unmatched records, and
(5) report paths. Every count with a scientifically meaningful denominator is
shown as `numerator / denominator (percentage)`; zero denominators are reported
as not applicable.

Concordance does not materialize the complete raw summary-statistics table and
extracted VCF together. It streams only the configured comparison columns once
into temporary Parquet staging tables while assigning immutable one-based
source-row IDs. The existing scientific comparison then runs sequentially for
each canonical chromosome value, including a separate invalid/null group.
Duplicate decisions and a single all-chromosome external EAF file use the same
partition key. A per-chromosome external EAF template is resolved with the VCF's
confirmed build and read only for the active chromosome, avoiding a combined
genome-wide reference table. Compact partition counts are combined globally,
while exact medians, Pearson and Spearman correlations are calculated by lazy
scans of the temporary matched records. Final reports retain deterministic
source-row ordering, and all staging files are removed on success or failure.
