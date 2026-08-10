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

At chromosome step 03, the supplied raw frequency reference remains an
ordinary `CHROM`, `POS`, `REF`, `ALT`, population-AF table. PostGWAS joins it
once to the study on chromosome and position, then uses vectorised Polars
expressions to classify the four orientations. It does not create or require a
four-times-expanded reference file. The selected result is recorded in
`strand_action`, with aligned reference ALT frequency in
`strand_reference_af`.

Complementing both alleles preserves beta, odds ratio, Z and true EAF. Changing
the effect allele negates beta and Z, reciprocates an odds ratio before its
later log conversion, and replaces true EAF with `1-EAF`; SE, p-value, INFO and
sample size remain unchanged. A confirmed MAF is never inverted because it is
not tied to the listed effect allele. If the second-level reference comparison
instead proves that a MAF-like column is true EAF, swapped rows are inverted at
that point. After chromosome processing, PostGWAS consolidates every completed
chromosome's `maf_reference_decision`. Only unanimous conclusive evidence
changes the final dataset frequency type. The manifest preserves the initial
MAF suspicion separately and records the final type, its source, and the
chromosome decision counts. Mixed, missing, or inconclusive evidence leaves the
final type unresolved. Palindromic candidates use study consensus and a true
study EAF versus the selected population AF; values in the configured
near-one-half band remain ambiguous. The defaults reject unmatched or
ambiguous rows with explicit provenance. After the final internal or external
EAF is known, `strand_af_difference` compares it with the aligned population
AF. Differences beyond `strand.af_tolerance` are warnings by default because
population frequencies are supporting evidence; YAML can instead reject or
fail them through `strand.af_discordance_action`. These decisions follow the
[GWAS Catalog summary-statistics harmonisation method](https://www.ebi.ac.uk/gwas/docs/methods/summary-statistics).

External EAF and INFO annotations use the same row-preserving direct/swapped
allele matcher after strand orientation. It normalizes the four join keys,
gives direct matches precedence, performs one Polars left join and preserves
the study row count and order. On a swapped match, EAF becomes
`1-AF`; INFO is allele-independent and remains unchanged. EAF retains its
separate palindromic-frequency rules, while INFO retains its independent range
and missing-value policies. A duplicated INFO reference either keeps the first
row under `info.deduplicate_reference: true` or fails before row multiplication.

### Effect scale before Z-based recovery

The study-wide effect decision is passed unchanged to every chromosome. At
chromosome step 05, an odds ratio is normalized through the existing
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
output root. Each dataset also receives its own run metadata and complete
dataset/chromosome logs. Success, permanent failure, configuration failure,
and user interruption all append an explicit final record; the screen remains
limited to progress, warnings, failures, and result locations.

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
under each dataset's `00_harmonised_sumstat/logs/adapters/gwas2vcf/` directory;
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

- `<dataset>_<build>_qc_assessment.tsv` contains raw and final metrics;
- `<dataset>_<build>_qc_filter_rules.tsv` contains raw-VCF rule and reason counts;
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
those values rather than maintained independently.

When the study-wide frequency distribution is MAF-like, chromosome step 03
checks it against the population column selected by `comparison_af.column` in
the 1000G table declared by `resource_layout.default_eaf`. Structural columns
and the tab delimiter come from `default_eaf_mapping`. Direct allele matches
use the ALT frequency; swapped matches use `1-AF`. The study column is accepted
when the listed effect allele is itself minor above the configured MAF-like
cutoff, or when input AF is at least as close to effect-allele-aligned reference
EAF. It is rejected when reference MAF is closer. Palindromic A/T and C/G
variants do not contribute to this decision.

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
complement, and complement-plus-swapped orientations. Indels require an exact
allele representation; the validator does not normalize them against a FASTA.
SNP and indel retention are therefore reported as separate strata. When every
configured SNP matching check passes and
`concordance_validation.indel_representation_action: warn`, balanced input-only
and VCF-only indel counts are reported as the maximum number of possible
normalization-related representation pairs within the same chromosome and
allele-length change; they do not by themselves fail the audit. This is a
count-based upper bound, not a record-level match: every such indel remains in
the input-only or VCF-only report, and unpaired indels still use the configured
failure thresholds. Set the action to `fail` to require
exact representation for all indels. Equal-length multi-base substitutions are
reported in a third `other_variants` stratum and are never treated as indels or
as possible normalization changes. This distinction is necessary because
[`bcftools norm -f`](https://samtools.github.io/bcftools/bcftools#norm)
left-aligns and normalizes indels before concordance.
It compares log-scale effects, effect-allele frequency, and Z score after
orienting all three to the VCF ALT allele. When the input has no Z column, Z is
calculated from effect/standard error, or from signed effect and the configured
p-value scale when standard error is unavailable. Palindromic comparisons and
all concordance tolerances are controlled by `concordance_validation` in the
harmonisation YAML. Input parsing and p-value interpretation reuse the canonical
`input` and `pvalue` policies rather than defining a second set of values.
Concordance reads the final reference-resolved frequency type from the run
manifest. It performs effect-allele-oriented EAF comparison after reference
confirmation of EAF, preventing `0.2` and an incorrectly inverted `0.8` from
appearing concordant after MAF folding. If the final type is unresolved, only
the allele-frequency comparison is skipped with a warning; variant matching,
effect, and Z-score checks still run.

Reports are written below
`<output>/<dataset>/harmonisation/00_harmonised_sumstat/qc_summary/concordance/`.
The status summary, value mismatches, input-only variants, and VCF-only variants
are always written; the full matched table is optional. A complete audit log is
written to the sibling `logs/` directory even when preflight or analysis fails.

Concordance does not materialize the complete raw summary-statistics table and
extracted VCF together. It streams only the configured comparison columns once
into temporary Parquet staging tables while assigning immutable one-based
source-row IDs. The existing scientific comparison then runs sequentially for
each canonical chromosome value, including a separate invalid/null group.
Duplicate decisions and external EAF data use the same partition key. Compact
partition counts are combined globally, while exact medians, Pearson and
Spearman correlations are calculated by lazy scans of the temporary matched
records. Final reports retain deterministic source-row ordering, and all
staging files are removed on success or failure.
