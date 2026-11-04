# Harmonisation configuration

## Canonical command

```console
postgwas harmonisation \
  --sample-sheet studies.csv \
  --run-config postgwas.yaml
```

Explicit CLI values override the user YAML, which overrides packaged YAML.
Argparse does not own configurable defaults.

For an internal frequency column confirmed as MAF, see the paired
[optional EAF replacement options](eaf_replacement.md).

## Sample sheet

The version-2 sample sheet contains dataset-specific paths, identifiers, and
source-column mappings. `resource_directory` and `output_directory` are run
settings and do not belong in dataset rows.

`external_eaf_file` and `external_info_file` accept either one existing table
containing all chromosomes or an explicit path template containing
`{chromosome}` and optionally `{build}`. A common filename prefix alone is not
resolved implicitly, because selecting a reference file by guesswork could use
the wrong chromosome or panel.

For a multi-chromosome study, a user table without `{chromosome}` is scanned
once at dataset step 08. PostGWAS projects only its configured chromosome,
position, allele, and EAF/INFO columns and writes one temporary Parquet
partition per observed chromosome before the process pool starts. EAF and INFO
columns in the same file and with the same structural mapping share that scan.
`external_reference_staging.batch_rows`, `compression`,
`compressed_suffixes`, and `atomic_output_suffix` are canonical YAML settings;
the temporary paths come from `output_layout.external_eaf_partition` and
`external_info_partition`.
Original values are not deduplicated, reoriented, or filtered during staging.
An input containing `{chromosome}`, and a run containing only one chromosome,
remain on the existing direct-read path.

`trait_type`, `effect_type`, `p_value_type`, and `delimiter` are inferable.
Omitted or missing values become `auto`; explicit `auto` remains supported.
Known values and aliases are case-insensitive and normalized to lowercase.
An unknown provided declaration fails during sample-sheet input validation,
with its row, field, original value and accepted choices. The CSV loader,
direct model and generated drafts share that validation. See the
[2026-10-31 behavior-change note](../../wiki/harmonisation/sample-sheet.md#behavior-change--2026-10-31)
when updating a sheet that previously relied on invalid values becoming `auto`.
Scientific inference and comparison with valid explicit values remain a later
preflight stage.

Explicit non-`auto` sample-sheet declarations have final precedence over the
generic YAML values for `effect.type`, `pvalue.type`, and `input.delimiter`.
PostGWAS still runs automatic inspection as an independent cross-check. When
the two disagree, `validation.declaration_mismatch_action` controls the result:
the default `warn` prints a bright-red terminal warning, records the mismatch,
and continues with the declaration; `fail` stops the dataset before chromosome
processing. Automatic inspection is evidence, not permission to silently
replace an explicit scientific declaration.

Summary-statistics parsing remains strict. A Python `MemoryError` during the
initial parse permits one retry with `low_memory=True`, using the same delimiter,
schema and error-handling options. Format, schema and other non-memory parser
errors propagate immediately rather than repeating the read. A failed retry
also propagates through the existing canonical failure logger. Lower-memory
parsing reduces memory pressure rather than repairing malformed data
([Polars reader documentation](https://docs.pola.rs/api/python/version/0.20/reference/api/polars.read_csv.html)).

### P-value representation

`pvalue.type` and the sample-sheet `p_value_type` accept exactly `auto`, `raw`,
or `neglog10`. These correspond to the two standard GWAS representations used
by the [GWAS Catalog summary-statistics format](https://www.ebi.ac.uk/gwas/docs/summary-statistics-format)
(`p_value` or `neg_log_10_p_value`) and the
[GWAS-VCF specification](https://github.com/MRCIEU/gwas-vcf-specification)
(`LP = -log10(p)`). Signed `log10(p)`, signed `ln(p)`, and positive `-ln(p)`
are not accepted; convert them to raw p or non-negative `-log10(p)` first.

Automatic detection runs once on the complete study before chromosome fan-out.
Neither supported representation can be negative. The dataset-level
`pvalue.max_negative_fraction` gate uses finite numeric cells as its denominator:
up to the packaged 0.1% threshold produces a prominent warning and rejects those
rows, while more than 0.1% stops before chromosome processing. A `neglog10`
decision then requires all three YAML-backed checks:
at least `pvalue.mlogp_detect_min_count` values and at least
`pvalue.mlogp_detect_proportion` of usable values must exceed
`pvalue.mlogp_detect_threshold`, and the study-wide median must be within
`pvalue.mlogp_median_tolerance` of `pvalue.mlogp_expected_median`. The canonical
YAML supplies the expected median as `0.3010299956639812` (`log10(2)`), because
a complete null-dominated GWAS has raw p median 0.5. A selected or top-hit-only
file may not retain this distribution and should declare `neglog10` explicitly
only when its source metadata proves that scale.

This design prevents one sentinel from choosing the scale: both the configured
count and fraction must agree, and even then the configured median check must
confirm the interpretation. Incompatible or insufficient evidence stops during
dataset step 07; no chromosome can receive a different decision.

After that decision, both supported representations produce one canonical raw
p-value column named by `pvalue.output_column` (packaged default `PVAL`). A raw
`0.05` remains `PVAL = 0.05`; a supplied `-log10(p)` value of approximately
`1.30103` is converted to the same `PVAL = 0.05`. Harmonisation does not create
a canonical `LP` companion. The GWAS-to-VCF adapter consumes raw `PVAL` and
writes the specification-defined `FORMAT/LP = -log10(p)` in the final VCF.
PostGWAS refuses to overwrite an unrelated input column that already has the
configured output name.

### Multiple INFO values and score type

A selected internal INFO column may contain a list such as
`0.95,0.80,0.60` in one variant row. `info.multi_value_delimiter` defines its
single-character separator and defaults to a comma. Dataset step 03 converts
that list to one scalar before duplicate validation.
`info.multi_value_aggregation` is the schema-validated YAML policy: `median`
is the shipped compatibility default, `mean` selects an arithmetic mean, and
`fail` requires the user to supply an already-combined scalar. Both numeric
choices are row-wise and unweighted.
Missing tokens are omitted; `info.multi_value_invalid_token_action` controls
whether a non-missing non-finite or non-numeric token is logged and ignored or
stops the dataset. The scalar is written to the collision-checked column named
by `info.multi_value_output_column`.

PostGWAS deliberately does not expose a sample-size-weighted mode for this
input representation. The sample sheet contains no ordered cohort-level N
vector that can be matched safely to an unlabelled ordered INFO list. A total
study N cannot reconstruct those weights. When equal cohort influence is not
appropriate, set the aggregation policy to `fail` and calculate a documented
scalar upstream from comparable cohort metrics and matching per-variant cohort
sample sizes. GWAS-SSF expects one numeric INFO value between 0 and 1, while
imputation-aware meta-analysis uses both study size and imputation information
when weighting study evidence ([GWAS Catalog
format](https://www.ebi.ac.uk/gwas/docs/summary-statistics-format),
[de Bakker et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC2782358/)).

The source column, selected policy, affected-row and invalid-token counts, and
working scalar column are written to the canonical log. The original text is
retained for duplicate evidence and rejected-row provenance. VCF metadata names
the original input column and records the unweighted aggregation. Because the
scalar exists before duplicate handling, the configured higher-INFO
tie-breaker is active for these studies rather than falling through to input
order.

`info.score_type` accepts `auto`, `standard_info`, or `mach_rsq` and defaults
to `auto`. The selected internal column, external column, or fixed value is
assessed once at dataset step 08 before chromosome fan-out. The resulting type
and evidence are stored in the run manifest and passed unchanged to every
worker; a chromosome cannot independently reinterpret the same study.

Automatic resolution uses finite numeric values from the complete selected
source. Values through `info.clip_tolerance` are standard-INFO evidence. Values
above that tolerance and no higher than `info.mach_rsq_max` are MaCH Rsq
evidence. `mach_rsq` is selected only when their fraction is at least
`info.auto_mach_rsq_fraction`; otherwise `standard_info` is selected and those
isolated high rows follow `info.out_of_range`. This prevents one 1.5 sentinel
in a large file from changing the treatment of every chromosome.

If a genuine MaCH Rsq source happens to contain no observed value above the
standard-INFO tolerance, the two scales are numerically indistinguishable.
Declare `info.score_type: mach_rsq` when source metadata establishes that scale
but the observed distribution contains no MaCH-specific evidence.

The shipped values are 1.05 for the standard-INFO tolerance, 2.0 for the MaCH
Rsq maximum, and 0.001 (0.1%) for the agreement fraction. These are canonical,
schema-validated YAML settings rather than literals in module code. On
`standard_info`, only values from `info.clip_max` through the tolerance are
treated as rounding overshoot and rescaled to `info.clip_max`. On `mach_rsq`,
values such as a chromosome-X Rsq of 1.5 are retained unchanged.

Values above `info.mach_rsq_max` are never evidence for either type. If their
finite-value fraction is strictly greater than `info.maximum_invalid_fraction`
(shipped as 0.001), the dataset stops before any worker starts and asks the
operator to verify the mapped column. At or below that guard the run continues
and chromosome step 11 applies `info.out_of_range`, whose shipped action is
`reject`. The former `clip` action is not available: a value such as 87 can no
longer be turned into a perfect score of 1.0. An explicit score type bypasses
automatic classification but not this gross-invalid guard.

This distinction follows the [GWAS Catalog standard-format INFO
range](https://www.ebi.ac.uk/gwas/docs/summary-statistics-format),
[PLINK's documented MaCH Rsq upper filter range](https://www.cog-genomics.org/plink/2.0/filter),
and [REGENIE's distinction between MaCH Rsq and IMPUTE
INFO](https://rgcgithub.github.io/regenie/options/). Post-merge
`modules.qc_summary.rules.info_max` remains an independently configured virtual
QC threshold. A stricter QC maximum can count a retained MaCH Rsq value as a QC
failure, but it does not remove that value from the delivered raw merged VCF.

### Chromosome scope and PLINK aliases

The packaged harmonisation policy retains autosomes 1–22 and chromosome X
only. The read-stage `chromosome.allowed` and post-split
`chromosome.allowed_after_split` defaults are identical, so Y, mitochondrial,
scaffold and alternative-contig rows cannot proceed to chromosome resource
preflight or VCF creation.

Chromosome labels are normalized before that filter. PLINK 23 becomes X, 24
becomes Y, 25 becomes X and 26 becomes MT; literal `XY`, `PAR1` and `PAR2` also
become X, while `M` becomes MT. Mapping PLINK 25/PAR to X is intentional:
[PLINK defines code 25 as XY/PAR and 26 as MT](https://www.cog-genomics.org/plink/2.0/filter),
and its [`--merge-x` and VCF export behavior](https://www.cog-genomics.org/plink/2.0/data)
represent PAR records on the X contig used by GRCh FASTA and VCF resources.
PostGWAS does not convert ordinary chromosome Y to X.

Coordinate QC separately records missing chromosomes and positions, Y and MT
exclusions, other excluded chromosome labels, and the number of PAR-labelled
rows normalized to X. Every removed row remains in the rejected-variant audit.

Use `--exclude-chr-x` to explicitly remove chromosome X:

```console
postgwas harmonisation \
  --sample-sheet /path/to/studies.csv \
  --resource-directory /path/to/resources \
  --output-directory /path/to/results \
  --exclude-chr-x
```

The equivalent run-configuration policy is:

```yaml
modules:
  harmonisation:
    policies:
      chromosome:
        exclude_x: true
```

The packaged value is `false`. Omitting the flag preserves a YAML value of
`true`; the flag takes no value and explicitly enables exclusion. The resolved
normalization and rename map determine which labels mean X. With packaged
settings, this includes `X`, `chrX`, `23`, `25`, `XY`, `PAR1`, and `PAR2`.
Exclusion takes precedence over configured chromosome lists without changing
the order or scope of the remaining chromosomes.

The reader preserves source-row IDs and original values, resolves chromosome
identity from separate or combined coordinates, then excludes X before
missingness, scientific numeric conversion, allele validation, duplicates,
build inference, strand consensus, frequency decisions, and chromosome work.
Malformed scientific cells confined to excluded X rows do not affect retained
rows. File readability, table structure, headers, and mappings must still be
valid. Even with positive `input.schema_inference_rows`, mapped positions and
scientific numeric columns are read as text before exclusion, then retained
rows undergo normal explicit conversion and QC. Structural delimiter,
row-width, and header errors still fail. Projected reference-table and streaming
reference reads likewise use text before exclusion, so X values cannot affect
reference validation, duplicate handling, or frequency screening. Excluded rows receive
`chromosome_excluded_by_user`, distinct from an
invalid or unsupported chromosome; raw-input counts include them once, while
scientific QC denominators use the retained scope. An X-only dataset fails
with an actionable message rather than producing a successful empty result.

X-specific references are omitted from startup and observed-chromosome
preflight, staging, workers, and optional `--validate` comparisons, including
when the configured reference inventory lists X. Shared resources needed for
retained chromosomes remain required. Changing `exclude_x` invalidates the
scientific checkpoint under the global resume policy; an identical completed
exclusion run can be reused after normal checkpoint validation. This option
belongs to the direct `postgwas harmonisation` command; harmonisation currently
has no pipeline adapter.

### Strand consensus and per-variant orientation

Dataset step 05 tests each candidate genome build with forward, swapped,
reverse-complement and reverse-complement-swapped allele matches. The same
coordinate joins provide the non-palindromic forward/reverse counts used at
dataset step 06, so the build references are not read a second time. A
consensus requires `strand.min_informative_variants` and the dominant fraction
configured by `strand.consensus_threshold`. Palindromic SNPs never contribute
because their allele letters cannot distinguish the competing strands.

Normal automatic build acceptance requires all four evidence checks: at
least `build.min_match_count` allele-compatible rows (default 100), at least
`build.min_reference_match_fraction` of the winning build's unique reference
markers on chromosomes present in the study (default 0.10), at least
`build.min_match_fraction` of coordinate-testable study rows matching the
winner (default 0.80), and at least `build.confidence_ratio` of all allele
matches supporting that build (default 0.90). A row is coordinate-testable
only when its exact chromosome and position occurs in at least one configured
build reference. Rows absent from both references are reported separately and
the winning build's input-wide match fraction remains in the log, but neither
is used as a coverage denominator because the build-check files are
representative marker panels rather than inventories of every possible study
variant.

Reference coverage is chromosome-scoped after the same chromosome
normalization is applied to study and reference labels. For example, a
chromosome-1 study with 964,783 variants matched 32,534 of 35,333 unique
relevant GRCh37 markers (92.08%), despite those matches representing only
3.37% of all study rows. The correct reference denominator accepts the strong
GRCh37 evidence instead of falsely requiring half of a dense GWAS file to
occur in a sparse marker panel. Each study row counts once per build and each
reference marker counts once for coverage despite exact duplicates. An
explicit `build.mode` remains authoritative and bypasses the automatic
evidence thresholds while still collecting strand evidence from the declared
build reference. Sampling informative SNPs for build inference is also used by
[MungeSumstats](https://al-murphy.github.io/MungeSumstats/reference/get_genome_build.html)
([Murphy et al., 2021](https://doi.org/10.1093/bioinformatics/btab665)); the
0.10 reference-coverage default is a configurable PostGWAS safety policy, not
a universal biological threshold.

#### Experimental fallback for low reference coverage

The canonical default `build.low_coverage_fallback_enabled: true` enables the
fallback after checking normal acceptance first. Set this policy to `false`
to require the normal coverage gate without fallback; an explicit value in an
existing run YAML takes precedence over the packaged default. The fallback can
waive **only reference coverage**, when it is the sole failed automatic gate.
The absolute match count, testable-row fraction, unique winner and winning-share checks must still
pass. Declared builds bypass automatic inference; normal automatic acceptance
is unchanged even if fallback-specific checks would fail.

The fallback compares exactly two references. Both panels must contain markers
on every coordinate-testable study chromosome; otherwise a zero competing count
could reflect missing reference scope. Missing scope is recorded and rejects
the fallback with a request for comparable panels. Entirely untestable
chromosomes do not enter this check and remain explicitly unassessed.

The default additional acceptance requirements are:

| Policy under `policies.build` | Default | Meaning |
|---|---:|---|
| `low_coverage_min_match_fraction` | 0.95 | Winner matches / coordinate-testable study rows, inclusive minimum. |
| `low_coverage_max_competitor_fraction` | 0.05 | Competitor matches / the same rows, inclusive maximum. |
| `low_coverage_min_unique_loci` | 10000 | Distinct winner-only compatible SNP loci, inclusive minimum. |
| `low_coverage_chromosome_min_competitor_loci` | 100 | Chromosome veto locus floor, combined with the fraction below. |
| `low_coverage_chromosome_min_competitor_fraction` | 0.95 | Chromosome veto competitor match fraction, inclusive minimum. |

Coordinate-testable rows include allele-incompatible hits and missing alleles.
The row fractions preserve duplicate study rows, with each row contributing at
most one match per build; a row matching both builds counts toward both
numerators. The distinct-locus floor counts canonical `(chromosome, position)`
once, requires different A/C/G/T single-base alleles compatible with the winner,
and excludes a locus if **any** study evidence there matches the competitor,
including an opposing allele pair on another row. Non-SNPs can supply ordinary
row matches but cannot increase either SNP-locus floor. SNP loci can be in LD;
distinct does not mean statistically independent. No study rows are removed.

Fallback is vetoed if any chromosome has at least 100 distinct competitor-only
SNP loci **and** competitor matches at least 95% of that chromosome's
coordinate-testable study rows. A sparse chromosome or an absent veto does not
certify consistency. Small mixed components can pass these gates.

This is **experimental dominant-build inference**, not proof that every row
uses one build. A 95% match fraction is not a posterior probability. Comparative
mapping is supported by [GWASLab's InferBuild documentation](https://cloufield.github.io/gwaslab/InferBuild/)
(which warns below 10,000 matches), [tidyGWAS](https://doi.org/10.1093/bioadv/vbaf262)
(which compares 10,000 sampled study rows), and the
[MungeSumstats upstream implementation](https://github.com/neurogenomics/MungeSumstats/blob/master/R/get_genome_build.R).
These sources do **not** validate the exact PostGWAS 95/5%, 10,000 winner-only
loci or chromosome 100/95% thresholds. They are configurable engineering
safety policies. The existing downstream strand checks remain in force.

For example, save the following as `build-fallback.yaml`, replacing the
filesystem placeholders with your paths. It uses the existing run-configuration
interface and exposes all fallback thresholds without changing ordinary gates:

```yaml
run:
  output_directory: /absolute/path/to/harmonisation-output
resources:
  root: /absolute/path/to/harmonisation-resources
modules:
  harmonisation:
    policies:
      build:
        mode: auto
        low_coverage_fallback_enabled: true
        low_coverage_min_match_fraction: 0.95
        low_coverage_max_competitor_fraction: 0.05
        low_coverage_min_unique_loci: 10000
        low_coverage_chromosome_min_competitor_loci: 100
        low_coverage_chromosome_min_competitor_fraction: 0.95
```

```bash
postgwas harmonisation --sample-sheet /absolute/path/to/studies.csv --run-config build-fallback.yaml
```

The returned `genome_build` evidence in the dataset manifest records
`acceptance_path` (`normal`, `declared`, `low_coverage_fallback` or `ambiguous`),
all failed normal gates, distinct SNP-locus counts, per-chromosome row counts and
fractions, scope omissions, vetoes and fallback gate outcomes. The canonical
log, study-wide report and HTML report retain this evidence, including on
failure. QC summaries label fallback acceptance as experimental. All resolved
policies participate in the existing scientific checkpoint identity; changing
any of them follows the global restart/resume policy.

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
Here “reference AF” means the configured population-frequency table; the FASTA
reference genome supplies bases but contains no population allele frequencies.

Chromosome step 03 has already converted an odds ratio and any declared raw
OR-scale SE to canonical log-odds BETA and SE. Complementing both alleles then
preserves BETA, Z and true EAF. Changing the effect allele negates BETA and Z
and replaces true EAF with `1-EAF`; log-scale SE, p-value, INFO and sample size
remain unchanged. Study AF acceptance now happens once, before chromosome processing.
The whole-study screen triggers a reference check only when more than
`eaf.maf_decision_cutoff` (default 95%) of usable study AF values are at or below
0.5. Otherwise the declared EAF contract is accepted without another frequency
concordance check. For a MAF-like distribution, PostGWAS first tests the study's
original effect-allele frequencies against `default_eaf` across the dataset.
Accepted EAF stays unchanged. Otherwise the same eligible pairs may support a
MAF interpretation; supported MAF requires the paired
[explicit EAF replacement source](eaf_replacement.md), while insufficient or
conflicting evidence stops the dataset. Chromosome workers reuse that decision
and apply ordinary allele swaps exactly once. The manifest keeps the initial
screen and final `eaf_reference_validation` evidence separately. The default
comparison itself neither filters rows nor changes their frequencies; an
explicitly selected replacement supplies the working EAF while preserving the
original input MAF and source snapshots. A palindromic A/T
or C/G SNP in `strand.mode: auto` first follows a strong full-study
non-palindromic strand consensus.
When that consensus is mixed, unresolved, or insufficient, a true internal EAF
tied to the study effect allele may resolve the row only if both study and
reference AF are outside the YAML-configured ambiguity interval, one candidate
is within `strand.palindromic_af_max_difference`, and it beats every competitor
by `strand.palindromic_af_min_error_margin`. A decisive frequency call that
conflicts with strong consensus is rejected as
`palindromic_frequency_conflict`; absent/unusable evidence is
`palindromic_orientation_unavailable`, excessive reference difference is
`palindromic_frequency_discordant`, and near-half or non-unique evidence is
`palindromic_ambiguous`. With `strand.mode: reference_aligned`, the operator
explicitly declares the selected reference strand, but PostGWAS still requires
each row to match direct or swapped reference REF/ALT; complement-only rows are
rejected. These frequency rules follow the principle
used by [TwoSampleMR](https://mrcieu.github.io/TwoSampleMR/articles/harmonise.html),
while the study-wide consensus follows the [GWAS Catalog orientation
method](https://ebispot.github.io/gwas-sumstats-harmoniser-documentation/Introduction/Orientation-of-palindromic-variants/).
The defaults reject unmatched or ambiguous rows with explicit provenance.
After the final internal or external EAF is known,
`strand_af_difference` compares it with the aligned population AF. A large
difference is supporting QC evidence, not a reason to reorient the row.
Non-palindromic discordance warns by default through
`strand.af_discordance_action`; palindromic discordance rejects by default
through `strand.palindromic_af_discordance_action`. These decisions follow the
[GWAS Catalog summary-statistics harmonisation method](https://www.ebi.ac.uk/gwas/docs/methods/summary-statistics).

External EAF and INFO annotations use the same row-preserving direct/swapped
allele matcher after strand orientation. Before alignment of an ordinary
sample-sheet external EAF source, its deduplicated raw frequency column is screened with
`eaf.maf_decision_cutoff`. This screen does not prove allele meaning and never
auto-converts MAF: a non-MAF-like result retains the declared external ALT/EAF
contract, while a MAF-like result must pass the independent comparison below.
An explicitly selected MAF replacement instead uses its declared ALT-frequency
semantics and bypasses this raw-frequency classifier; numerical, coverage and
allele-matching safeguards still apply. `default_eaf` remains the sole source
of comparison and orientation evidence.
The matcher normalizes the four join keys,
gives direct matches precedence, performs one Polars left join and preserves
the study row count and order. On an ordinary swapped match, EAF becomes
`1-AF`; INFO is allele-independent and remains unchanged. For a palindromic
external EAF, the already-resolved `strand_action` is required. The external
mapping declares its REF and ALT/effect-allele columns, so the direct/swapped
join supplies ALT-aligned AF as listed or as `1-AF`; the panel is not allowed
to re-decide study strand. Missing upstream orientation is rejected as
`palindromic_orientation_unavailable`. This external comparison never changes
`strand_action`, BETA, Z, or the study effect direction, because an independent
external frequency is not provenance for the study beta. Chromosome QC records
EAF provenance as `study_supplied` or `reference_imputed`. INFO retains its
independent range and missing-value policies. The
strand, external EAF, and external INFO references use the single
`external_reference` duplicate policy. The defaults keep one row only when the
normalized chromosome, position, alleles and selected numeric annotation are
identical (`exact_duplicate_action: keep_one`). If values differ—including a
missing value beside a finite value—all reference rows for that allele key are
discarded (`non_identical_duplicate_action: discard_all`), so row order never
chooses scientific data. The affected study row is then unmatched and follows
the existing EAF or INFO missing-value policy. Duplicate-group and removed-row
counts are recorded in the chromosome log.

Frequency descriptions reuse earlier scalar statistics only within the current
validation/harmonisation call. A cheap cloned Series retains the selected
column's previous state. Reuse requires the same selected name, dtype, ordered
values/nulls, chunk lengths and effective `eaf.clip_tolerance`. Changed rows,
allele-frequency flips, clipping/nulling, source-column changes or rechunking
therefore cause the existing summary helper to recalculate. Zero-valued
frequencies are conservatively recalculated because numeric equality does not
distinguish positive and negative zero in descriptive metadata. Independent
before/after dictionaries, all scientific checks, decisions, warnings,
rejection provenance and report labels remain intact. No summaries are cached
across chromosomes, datasets or invocations, and no new tolerance is introduced.
This uses the public Polars [Series clone](https://docs.pola.rs/api/python/stable/reference/series/api/polars.Series.clone.html),
[equality](https://docs.pola.rs/api/python/stable/reference/series/api/polars.Series.equals.html)
and [chunk-length](https://docs.pola.rs/api/python/stable/reference/series/api/polars.Series.chunk_lengths.html)
contracts; a full-column hash or a second statistics implementation is not used.

Standalone frequency validation uses the resolved `eaf.maf_decision_cutoff`,
`eaf.invalid_fraction_cutoff` and `eaf.missing_fraction_cutoff` when threshold
arguments are omitted or `None`. Explicit helper arguments, including zero,
retain precedence. The chromosome workflow delegates to that same resolution.

### Effect scale before Z-based recovery

The study-wide effect decision is passed unchanged to every chromosome. At
chromosome step 03, an odds ratio is normalized through the existing
`harmonise_effect_estimates()` implementation, producing the canonical
log-odds beta and applying the configured non-positive-OR and SE-scale rules.
Only then does step 06 derive a missing standard error. It first requires the
harmonised BETA and signed Z to have the same non-zero sign, then calculates
the positive value `SE = abs(BETA/Z)`. Thus an OR uses `abs(ln(OR)/Z)`, never
`OR/Z`, and the canonical column mapping prevents a second logarithm later in
the workflow. Dataset step 07 first measures this consistency once across the
complete study. Its denominator contains only rows that need SE recovery and
have a finite valid effect plus a finite non-zero signed Z; missing values and
populated study SE cells cannot dilute it. With the packaged
`effect_from_z.max_beta_z_sign_mismatch_fraction: 0.01`, a discordant fraction
equal to 1% continues, while a fraction strictly above 1% stops before
chromosome partitioning or worker launch. This threshold is configurable and
schema validated rather than embedded in the algorithm.

At or below that limit, the default
`effect_from_z.beta_z_sign_mismatch: reject` records opposite signs—or zero
BETA with non-zero Z—as `beta_z_sign_discordant` and removes those variants per
chromosome. The screen and canonical log report both the complete-study
fraction and chromosome rejection counts. The stricter `fail` action stops on
the first non-zero complete-study count. Populated study SE values are not
changed by this policy. If the input has no effect column, step 05 is skipped
and step 06 creates beta directly from Z.
This follows the standard logistic-regression convention that the reported
coefficient is the log odds ratio and its Z statistic is coefficient divided
by standard error ([Stata logistic-regression interpretation](https://www.stata.com/links/stata-basics/logistic-regression-1-introduction/)).

When neither BETA nor SE was supplied, step 06 reconstructs both from Z, EAF,
and Neff. The default `effect_from_z.method: metal_large_n` uses
`BETA=Z/sqrt(2*EAF*(1-EAF)*Neff)` and
`SE=1/sqrt(2*EAF*(1-EAF)*Neff)`. It is the large-N, standardized-trait
approximation intended for sample-size-weighted Z results such as METAL output.
Relevant primary references are [METAL (Willer et al. 2010, PMID
20616382)](https://pubmed.ncbi.nlm.nih.gov/20616382/) and the educational-
attainment meta-analyses of [Rietveld et al. 2013 (PMID
23722424)](https://pubmed.ncbi.nlm.nih.gov/23722424/), [Okbay et al. 2016 (PMID
27225129)](https://pubmed.ncbi.nlm.nih.gov/27225129/), and [Lee et al. 2018
(PMID 30038396)](https://pubmed.ncbi.nlm.nih.gov/30038396/). These references
provide the sample-size-weighted meta-analysis context; PostGWAS records the
large-N conversion as an approximation rather than claiming to recover a
study's fitted coefficient.

The alternative `effect_from_z.method: zhu_2016` uses the finite-sample
adjustment `Neff+Z^2` in both denominators, following [Zhu et al. 2016 (PMID
27019110)](https://pubmed.ncbi.nlm.nih.gov/27019110/). With the default
`effect_from_z.phenotype_standard_deviation: null`, both methods assume a unit
phenotype standard deviation and therefore return standardized additive
effects. Supplying a positive phenotype standard deviation multiplies BETA and
SE by that value, preserving Z while putting a quantitative-trait approximation
in those phenotype units. This does not recover a logistic log odds ratio or a
liability-scale effect for a binary trait.

The diploid variance `2*EAF*(1-EAF)` is not a generally valid chromosome-X
model. Male non-PAR dosages may be coded `0/1` or `0/2`, PAR and non-PAR loci
have different ploidy, and a pooled EAF and Neff do not recover sex-specific
sample composition or allele frequencies. Therefore
`effect_from_z.x_chromosome_z_only_action` defaults to `fail`: when X is
present and the study supplies Z but neither BETA nor SE, the dataset stops
before genome-build/strand inference, resource preflight, chromosome
partitioning, or worker launch. No X effect is calculated. Supplying BETA or
SE remains allowed because sign-consistent `SE=abs(BETA/Z)` and `BETA=Z*SE`
do not use an autosomal genotype-variance model.

Two alternatives are explicit. To produce an autosome-only output, set both
`chromosome.allowed` and `chromosome.allowed_after_split` to lists containing
chromosomes 1 through 22 only; the normal input-rejection audit then records
excluded X rows. To reproduce the previous
numeric approximation, set
`effect_from_z.x_chromosome_z_only_action: allow_autosomal_assumption`. That
override produces a prominent warning and is recorded in chromosome QC and
the merged VCF effect/SE provenance; it is not described as an X-specific or
study-scale reconstruction. PostGWAS does not make an automatic PAR exception
because its canonical input convention maps `PAR1`, `PAR2`, `XY`, and PLINK 25
to the X contig. This policy follows the model-dependence described by
[Clayton 2008](https://pubmed.ncbi.nlm.nih.gov/18441336/), the alternative male
dosage models exposed by [PLINK 2](https://www.cog-genomics.org/plink/2.0/assoc),
and sex-chromosome summary-statistics QC guidance from
[Winkler et al. 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4083217/).

PostGWAS never divides unless the selected denominator is finite and positive.
Exact endpoint frequencies cannot be used even when `eaf.degenerate: keep`
permits them beside supplied effect statistics. Before calculating, PostGWAS
also requires the configurable proxy `2*EAF*(1-EAF)*Neff` to meet
`effect_from_z.minimum_effective_variance` (default `1.0`). The default action
rejects low-information rows with a distinct reason; `fail` stops the
chromosome. This proxy is not an exact minor-allele count because Neff may
differ from total N and EAF may come from a reference panel.

The [GWAS Catalog summary-statistics format](https://www.ebi.ac.uk/gwas/docs/summary-statistics-format)
defines EAF as a probability between 0 and 1. Non-finite EAF and finite values
outside `[0,1]` are rejected by default as
`eaf_out_of_range`. Clipping remains available only as an explicit policy; it
is no longer applied after `null`, `reject`, or `fail`. Exact `0` and `1` are
handled separately by `eaf.degenerate` so malformed inputs and probability
endpoints retain different provenance.

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

Chromosome step 09 only produces Z: it retains a supplied value or calculates
`Z = BETA / SE` where BETA and SE are finite and SE exceeds the positive,
configurable division floor. An unusable calculation input leaves Z empty; no
validation policy is applied and no row is removed in this step.

Chromosome step 10 owns the six basic checks and runs them exactly once: SE
missing, SE non-finite, SE at or below the division floor, BETA missing or
non-finite, BETA exactly zero, and Z missing or non-finite. The SE masks are
mutually exclusive and finiteness is tested before the denominator threshold,
so positive infinity cannot pass and one row is not counted under multiple SE
reasons. There is no DataFrame-identity certificate or conditional second pass.

The same step then validates the standard Wald relationship `Z = BETA / SE`
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

The resolved integer Neff is exported through the required `neff_col` adapter
mapping and written unchanged to `FORMAT/NEF`. The adapter does not recalculate
or round it again. The existing sample-size formula and PostGWAS rounding are
unchanged; cases, controls and total N retain their separate fields. The
selected Neff column and its statistics are included in the adapter export
summary. Missing or invalid adapter Neff values fail rather than being replaced
with a recalculation. This avoids a second rounding rule such as Python's
[half-to-even `round`](https://docs.python.org/3/library/functions.html#round).

Initial command preflight validates all sheet declarations, required source
pairs, path-template syntax, duplicate identities and shared setup without
opening the study tables. Each dataset's input availability, header, values
and active external EAF/INFO resources are checked when its turn begins.
The dataset then completes the full harmonisation pipeline before the next
study is read. A dataset-specific failure is recorded while independent later
datasets continue. Invalid configuration or sheet declarations, shared resource
failures, mandatory reporting failures and user interruption stop the batch.
After all independent datasets have been attempted, any failed or partial
dataset makes the command exit nonzero. Successful dataset outputs remain
available, but the incomplete batch cannot receive a completed checkpoint.
The direct sample-sheet loader retains its strict default file-metadata checks.

## Execution defaults

Packaged automatic rules reserve one logical CPU and ten percent of physical
memory. The actual resolved values are shown in `--help` and written to
`resolved_config.yaml`.

The shared `--threads` and `--memory-gb` options are identical for standalone
modules and pipeline mode. Their implementation lives in one common compute
CLI component. No deprecated harmonisation interface is available.

For chromosome processing, `--threads` is a total CPU budget rather than a
worker count. Adaptive scheduling is the default. It allows up to
`floor(threads / execution.min_threads_per_chromosome)` workers (default minimum
2), subject to chromosome count and memory. The allocated threads never exceed
`execution.threads_per_chromosome` (default maximum 5); a one-thread CPU budget
uses one worker with one thread. Polars and bcftools receive the same allocation.
PostGWAS places
`POLARS_MAX_THREADS` in the inherited environment before the spawn pool starts
and verifies the resulting pool in the child initializer. This order matters
because [Polars documents the pool as immutable after process start](https://docs.pola.rs/api/python/stable/reference/api/polars.thread_pool_size.html).

Dataset step 08 writes the validated study one chromosome at a time instead of
materialising every chromosome partition together. A compact lookup groups the
working row positions once and matches immutable source positions once by the
stable one-to-one input-row ID. It retains native position-index lists rather
than copies of every chromosome's scientific columns. Each chromosome's working
rows are selected in their existing order; its source rows are selected in
stable source-ID order. IDs are never interpreted as physical row offsets.
Actual index and ID order is checked before sorting, so already ordered data
avoids redundant sorts without relying on a backend's default join order.
Verified consecutive physical positions use a
[slice](https://docs.pola.rs/api/python/version/0.20/reference/dataframe/api/polars.DataFrame.slice.html);
ascending scattered positions use a temporary Boolean mask, and differently
ordered positions use indexed gathering. The mask uses
[NumPy](https://numpy.org/doc/stable/reference/generated/numpy.zeros.html)
and Polars' order-preserving
[filter](https://docs.pola.rs/api/python/version/0.20/reference/dataframe/api/polars.DataFrame.filter.html).
Native bitmap filtering still materialises each chromosome, while chromosome
assignment and source-ID matching are computed once. Only one chromosome's
subsets and a temporary mask are materialised at a time.
This uses Polars [row indices](https://docs.pola.rs/api/python/version/0.20/reference/dataframe/api/polars.DataFrame.with_row_index.html)
and [ordered grouping](https://docs.pola.rs/api/python/version/0.20/reference/dataframe/api/polars.DataFrame.group_by.html)
without changing row membership or scientific values. Both subsets are atomically
written as typed Parquet before the next chromosome is selected. The worker therefore
recovers the exact validated schema without a text round trip or a second type
inference. `input.chromosome_partition_compression` controls only this internal
storage compression (default `zstd`); it does not alter values, row order,
chromosome assignment, or reject provenance. Chromosome workers still start in
parallel after partitioning and resource staging have completed.

After partitioning, adaptive mode reads only Parquet metadata and reference
file sizes and reuses validated chromosome row counts. It estimates each
worker's memory as:

`memory_safety_factor * (worker_memory_floor_gb + working_data_gb + reference_gb + sort_buffers_gb)`

Working data is the larger of `rows * memory_bytes_per_variant` and the
study/source Parquet metadata size times `study_memory_multiplier`. References
named by `memory_reference_keys` use `reference_memory_multiplier` times their
on-disk sizes, or uncompressed Parquet metadata sizes. Paths reused for two
roles are counted once per chromosome. Metadata is cached within planning.
The defaults are 1 GiB base, 4096 bytes per variant, 8-fold study/reference
multipliers and a 1.5 safety factor. These are configurable planning estimates,
not measured peaks or guarantees. In particular, encoded Parquet sizes and
compressed text sizes do not equal resident memory.

Workers share the resolved memory budget after reserving the larger of parent
RSS and `memory_headroom_fraction` (default 0.10). Admission alternates the
largest and smallest remaining jobs that fit. Only admitted jobs are submitted;
when one completes, its reservation is released and the next fitting job starts.
A chromosome whose estimate cannot fit even alone fails planning with an
actionable message. OOM failures still stop the dataset without identical retries.
Scientific transformations and output merge order are unchanged.

Set `modules.harmonisation.policies.execution.scheduling_mode: fixed` in a run
configuration to use the previous fixed 20 GiB reservation and 5-thread allocation.
Adaptive controls live under that same execution policy block. Direct library
callers selecting adaptive mode must supply `total_cpu_budget` and
`memory_budget_gb`; public direct and pipeline commands resolve these globally.
The dataset log and run manifest retain estimates, parent headroom and CPU limits;
the log records admitted chromosomes and their active reservations. Benchmark
peak process-tree RAM and elapsed time on representative datasets before tuning
the estimate factors down. Pool allocators may retain memory between jobs, and
unusually wide tables or highly duplicated panels can exceed these estimates.

Implementation references: [Python futures](https://docs.python.org/3/library/concurrent.futures.html)
and [Arrow Parquet metadata](https://arrow.apache.org/docs/python/generated/pyarrow.parquet.RowGroupMetaData.html).

## Runtime prerequisites

### Build-liftover allele handling

`vcf.liftover_swap` has two values. The shipped `exclude` default counts and
removes successfully harmonised records marked `SWAP=1` or `SWAP=-1`; `keep`
retains them. Both modes use the same mandatory allele-aware correction. There
is no separate `update_tags` mode and no free-form plugin-tag argument.

Before liftover, PostGWAS reads the annotated VCF header. It validates the
canonical study fields `INFO/AF`, `FORMAT/AF`, `FORMAT/ES`, and `FORMAT/EZ`,
then derives population-frequency INFO fields from
`vcf_processing.external_frequency_columns`. Every field passed as an AF or
signed-effect tag must be declared with `Number=A` and `Type=Float`. The
resolved lists are written to the canonical log. The packaged configuration
does not include `FORMAT/AP1`, `FORMAT/AP2`, or `FORMAT/ED`.

For an allele swap, the plugin reorders ALT-frequency values and reverses the
sign of ES and EZ so that they continue to describe the resulting ALT allele.
This follows the allele-aware behavior described by the
[bcftools/liftover publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC10832354/)
and the [pinned upstream plugin](https://github.com/freeseek/score).

Target-build records are normalized and split, globally sorted, and only then
deduplicated with `bcftools norm -d exact`. The log records plugin input,
plugin rejections, swapped and reference-added records, policy exclusions,
normalization changes, duplicates removed, and final output separately.

Harmonisation requires `bcftools`, `tabix`, and the bcftools `liftover` plugin.
The project installs bcftools/HTSlib 1.23.1 from the environment and compiles
the plugin against matching bcftools source. Both the release archive and the
immutable plugin source revision are SHA-256 verified. Check an installation
with:

```console
bcftools --version
bcftools plugin -l | grep '^liftover$'
```

The recommended `tools/setup/install_postgwas.sh` installer builds the pinned
plugin inside the environment automatically. For a manual installation, first
activate the environment and run `tools/setup/install_bcftools_liftover.sh`.
Optional first and second arguments override the bcftools version and plugin
source revision; arguments three and four must then provide the corresponding
expected SHA-256 values.

The command verifies these requirements before reading a summary-statistics
file. A missing program or plugin therefore produces an actionable preflight
error in `run_metadata/postgwas.log`, rather than failing after chromosomes
have already been processed.

## Field completeness and recovery reporting

Dataset step 03 applies `columns.mandatory` as a row-level scientific rule. A
plain field is required independently; a nested group such as `[beta, zscore]`
requires at least one member. EAF, SE, INFO, SNP, and effective N remain absent
from the packaged mandatory list because the canonical `field_lifecycle` table
names their later source or derivation step. A row rejected at this gate uses
the explicit `missing_read_mandatory_value` reject reason.

The canonical `output_layout.field_completeness` path writes one structured TSV
per dataset. It records, without hard-coded field thresholds, the resolved study
column, post-parse missing count and fraction, read requirement, selected input
alternative, recovery step, final-check requirement, final action, and
post-recovery missing count. BETA and Z also receive an explicit alternative-
group row showing how many variants lacked both. The same structure is retained
in the run manifest and summarized in the canonical log and screen report.

Input missingness is measured after explicit chromosome exclusion and
configured missing-token parsing but before
numeric normalization; later unparseable numeric values are reported by the
normalization and field-specific validation steps. Final missing counts are
measured after recovery at chromosome step 13. With the packaged `reject`
action, they are attributed sequentially in `final_check.require` order, so a
row missing several final fields is counted under the first configured rule
only. With `keep`, the table instead labels the counts as independent and
potentially overlapping; `fail` is labelled as fail-fast. Chromosome coverage
is recorded so a partial run is not presented as a complete genome-wide
assessment. An ordinary sample-sheet external EAF source is selected only when
no internal EAF column is selected; it does not silently fill individual nulls
in an internal study column. The separate
[explicit EAF replacement](eaf_replacement.md) can replace the working frequency
source only after whole-study reference evidence supports MAF. SE has no
external-file source. When no SE column is selected, it is
derived from usable Z or from BETA and P under the configured scientific
policies. When a selected study SE column is only partly populated, the
pipeline preserves every populated value, first recovers eligible null cells
exactly from BETA and Z, and then uses the packaged
`pvalue.derive_partial_missing_se: true` BETA/P fallback for null cells still
unresolved. The p-based inverse-normal calculation uses preserved `ln(P)` and
`scipy.special.ndtri_exp`, so a positive raw token such as `1e-400` or a
negative-log10 value of `400` is not clipped before SE recovery. The derived
counts are recorded in chromosome QC, and the recovery hierarchy is recorded
in VCF provenance. Setting the policy to `false` disables only the p-based
fallback and leaves remaining null cells for `validation.se_invalid`; literal
input `p=0` continues to follow `pvalue.zero_missing_se`.

When supplied SE is complete, or ordinary partial recovery is disabled and
no missing-SE zero-P case needs attention, the SE step skips constructing
log-P provenance. Input validation, numeric normalization and SE diagnostics
still run. Missing-cell recovery and explicit zero-P handling retain their
existing log-domain calculations, source provenance and policies.

The standalone BETA/P-to-SE helper also uses `pvalue.clip_high` when the range
action clips an invalid raw `P>1`. It does not replace a valid `P=1` with an
independent approximation. In the default two-sided test, `P=1` gives `Z=0`,
so a finite SE cannot be recovered and remains missing for the existing gate.
Clipping provenance and the configured clipped-SE validation remain in force.

## Logs and run metadata

Every run writes `run_metadata/postgwas.log` and `resolved_config.yaml` at the
output root. It also writes the configured `runtime.run_summary_file` (packaged
default `harmonisation_run_summary.csv`) with one row per selected dataset.
Its columns follow pipeline order. This run-level CSV records final status,
failure reason, recorded failure stage and dataset log path;
input-QC-ready SNP and indel/other counts and `excluded_chromosome_variants`;
resolved and detected build/effect/P
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
The final terminal summary names failed, partial and unreached datasets with
their recorded reasons and available logs. The stage names describe the
observed orchestration boundary; they do not infer a scientific failure stage
from exception wording. No configuration key or scientific default is added.

The configured `output_layout.html_report` writes one standalone HTML report
for each selected dataset (packaged pattern
`{dataset_id}_harmonisation_report.html`). The configured
`runtime.run_html_report` writes the combined report (packaged filename
`harmonisation_run_report.html`) beside the run-summary CSV. Both are generated
from the existing run-summary records and manifests. No summary-statistics or
VCF file is reopened and no new scientific metric is calculated. HTML output
is atomic and required in the same way as the run-summary CSV.

Normal terminal output is always copied incrementally to the configured
`output_layout.screen_report`. With the packaged layout, each dataset receives
`<output>/<dataset_id>/harmonisation/<dataset_id>_screen_report.txt`. Shared
run preparation and final summary blocks are copied to every selected dataset;
dataset processing blocks are written only to the dataset they describe.
The shared `logging.show_screen` setting defaults to `true`. Use
`--hide-screen` to suppress normal terminal output while continuing to write
both these dataset reports and the shared `logging.screen_log_file`. Set
`logging.show_screen: true` in YAML to enable display when a run configuration
has disabled it.

Each dataset also receives its own run metadata and complete dataset/chromosome
logs. Success, permanent failure, configuration failure, and user interruption
all append an explicit final record; the screen remains limited to progress,
warnings, failures, and result locations.

Dataset run manifests retain their existing JSON format and are published by
the shared atomic report writer. Serialization and writing finish in a
temporary file in the destination directory before replacement; failure before
publication preserves the last complete manifest and removes the temporary.

Worker-failure displays read logs incrementally and retain only the requested
tail (the existing default is 400 logical lines), with the exact omitted-line
count. Complete diagnostics remain in the source logs; a long diagnostic line
is retained intact, so tail memory includes the retained lines and the longest
physical input line. The combined dataset/chromosome transcript uses streaming
copy and the shared transactional publication helper. Source/output aliases
are refused, and a failed copy preserves any previous complete transcript.
The snapshot is written after the final post-merge step exits, including its
cleanup messages and final `STATUS` record. Source logs remain authoritative
for early failures and any later outer-finalization errors.

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

Duplicate handling has two passes. The first runs during dataset-level input
validation before chromosome partitioning. `duplicates.key` identifies groups;
its default is chromosome, position, effect allele, and other allele. Allele
letter case is ignored. When position and both alleles participate, PostGWAS
first removes shared trailing and then leading VCF padding in private key
columns, adjusting only the private key position for a removed leading base,
and normalizes the result to `min(EA,OA), max(EA,OA)`. Original study columns are
unchanged. Therefore A/G and G/A at one coordinate identify the same physical
allele pair, padded `ATT/AT` and minimal `AT/A` identify the same indel, while
A/G and A/T remain distinct variants.

This key-only rule follows the sequence-allele padding contract in the
[VCF specification](https://samtools.github.io/hts-specs/VCFv4.5.pdf) and the
representation normalization performed by
[`bcftools norm`](https://samtools.github.io/bcftools/bcftools#norm). It is
deliberately narrower than reference-backed left-alignment.

A normalized group containing both EA/OA orders is classified as
`swapped_orientation`, counted separately, and removed in full with rejection
reason `swapped_orientation_duplicate`. PostGWAS does not reconcile such rows
at this early stage: BETA or odds ratio, signed Z, and EAF are tied to the
reported effect allele, while effect type, frequency semantics, strand, and
reference orientation have not yet been resolved. Even numerically equal
swapped rows are therefore unsafe to collapse. `duplicates.conflicting_action`
controls whether processing continues after removal (default) or the dataset
stops after both reports are written.

For every same-orientation duplicate group, PostGWAS compares the non-empty
columns selected by `duplicates.consistency_fields`. If those scientific values
agree, one row is selected using the explicit `duplicates.selection_order`:
completeness across configured quality fields, larger per-variant effective
sample size, higher INFO, and finally the earlier physical input row.
When an internal INFO cell contains delimited values, its configured
unweighted scalar aggregation is completed first, so this INFO criterion uses
the actual scalar rather than treating the raw text as unavailable.
Case-control effective sample size uses `4/(1/Ncase + 1/Ncontrol)`; a
quantitative-trait N is used directly. The formula is the usual balanced
case-control equivalent sample size
([Mallard et al.](https://onlinelibrary.wiley.com/doi/full/10.1002/gepi.22609)).
The input-order tie-breaker makes the result reproducible without selecting on
statistical significance. Only duplicate rows are sorted for ranking; the
complete GWAS frame is not globally sorted and its input order is restored.

If a group contains different non-empty association values, every row in that
group is rejected as `conflicting_duplicate`; `duplicates.conflicting_action`
can instead stop the dataset after the reports are written. The duplicate TSV
contains every member of each group plus its `consistent`, `conflicting`, or
`swapped_orientation` classification, conflicting fields, and
`duplicate_action`. The row retained for analysis is labelled `kept`; other
consistent rows are labelled `removed`, while unsafe rows carry the configured
conflict action. The canonical log reports swapped-orientation group and row
counts separately. The rejected-variant file records every swapped row under
`swapped_orientation_duplicate` and every same-orientation conflict under
`conflicting_duplicate`. `duplicate_input_row` is the one-based physical source
row used by concordance to identify that exact row; concordance does not repeat
the duplicate ranking.

Rejected-variant provenance uses the same stable one-based source-row ID. The
input reader captures an immutable snapshot immediately after parsing and before
coordinate, allele, frequency, effect, or standard-error transformations. Only
the snapshot for the active chromosome is loaded when rejected rows are written,
so the rejected file contains the original parsed study values even when the
working row was subsequently swapped, inverted, rescaled, or converted. The
internal chromosome snapshots are compressed Parquet intermediates and are
removed during successful final cleanup.

This policy deliberately does not keep the smallest p-value from conflicting
records, because selecting on significance can bias the retained result. The
[GWAS Catalog format](https://www.ebi.ac.uk/gwas/docs/summary-statistics-format)
permits either allele order but defines effect and frequency fields relative to
the effect allele, supporting an unordered physical key plus conservative
handling of mixed orientations. MungeSumstats separately documents duplicate
and allele-flip checks; its upstream
[`check_dup_row`](https://github.com/Al-Murphy/MungeSumstats/blob/master/R/check_dup_row.R)
uses an ordered A1/A2 key and keeps one duplicated key, whereas PostGWAS records
the unsafe orientation explicitly and does not choose a survivor before
orientation. Additional conventions are described by the
[GWAS-VCF publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC7805039/), the
[GWAS Catalog summary-statistics format](https://www.ebi.ac.uk/gwas/docs/summary-statistics-format),
and the
[MungeSumstats duplicate checks](https://www.bioconductor.org/packages/release/bioc/vignettes/MungeSumstats/inst/doc/MungeSumstats.html).

The second pass runs after reference orientation and final effect-frequency
alignment are complete within chromosome step 04. At that boundary the effect
allele is ALT, the other allele is REF, and BETA, Z, OR, and study EAF have been
transformed where required. Opposite-strand input rows such as A/G and T/C can
therefore converge on one reference-defined REF/ALT identity. The same
configured consistency fields and deterministic quality ranking then apply.
Finite aligned numeric values use the configured
`duplicates.post_orientation_relative_tolerance` to ignore only floating-point
roundoff introduced by transformations such as `1-EAF`. Compatible groups
retain one row and record the others as
`duplicate_variant`; disagreement between already aligned statistics rejects
the complete group as `conflicting_duplicate`. Each chromosome writes its full
post-orientation group evidence to the configured
`output_layout.post_orientation_duplicates` table, while the rejected-variant
output preserves the original input rows and the step-04 reason.

This second pass does not reverse-complement allele strings to construct an
indel key. Minimal key trimming is not FASTA-backed left-alignment and therefore
does not collapse repeat-shifted indels. Indels are otherwise compared only
after the supplied reference has matched their direct or swapped REF/ALT
representation; later bcftools normalization continues to left-align them. With
the canonical duplicate key, different ALT alleles at the same coordinate
remain different variants because both REF and ALT stay in the key.

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
After every required merged group and the final QC outputs validate, an `OK`
run removes the raw merged gwas2vcf VCF and its exact `.tbi`/`.csi` index by
default. Set `vcf.keep_gwas2vcf_intermediate: true` or pass
`--keep_gwas2vcf_intermediate` to retain and return it. This storage policy does
not change either annotated build-specific VCF or any variant/statistical value.
Cleanup runs only after post-merge status is `OK`, so earlier failed or partial
finalization can leave the raw merge available for diagnosis.
Validated inputs are ordered by `chromosome.allowed_after_split` before
concatenation, and that configured order must agree with the reference-contig
order. Merge groups use the resolved source and target builds plus configured
output-layout paths; build-like text in a dataset ID is never used to classify
a file. This enforces the [VCF 4.5 ordering contract](https://samtools.github.io/hts-specs/VCFv4.5.pdf),
which requires each chromosome to form one contiguous block with increasing
positions, and the [bcftools concat contract](https://samtools.github.io/bcftools/bcftools#concat),
which requires chromosome files to be supplied in output sort order.

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

The QC module then runs its shared `bcftools query` data pass against the
input-build merged VCF. It infers that file's coordinate build from its exact
`##genome_build=<build>` header declaration and extracts only the configured
fields into a temporary TSV. Two streaming Polars aggregation passes over that
table perform four related calculations without rereading the VCF. The first
pass collects the raw and QC-passed effective-sample-size mean and sample
standard deviation plus the configured raw reference quantile; the second uses
those scalars for all remaining metrics:

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
`modules.qc_summary.rules.sample_size_outlier_standard_deviations` standard
deviations. The default is
5, matching
[MungeSumstats `N_std`](https://www.bioconductor.org/packages/release/bioc/manuals/MungeSumstats/man/MungeSumstats.pdf).
This is descriptive QC: unusually large Neff values are reported but are not
silently removed and do not set the statistical-quality warning status. The
upper-tail count is retained only as an informational distribution diagnostic.

PostGWAS also reports the scientifically relevant lower tail. It calculates a
single reference, using linear quantile interpolation, from usable Neff in the
raw merged VCF and multiplies it by
`modules.qc_summary.rules.sample_size_minimum_fraction_of_reference`. The
packaged `sample_size_reference_quantile: 0.9` and two-thirds minimum fraction
match the default LDSC convention of the pandas-default 90th percentile divided
by 1.5. The same resulting threshold is applied to the raw and virtual
QC-passed stages, and each stage reports the count and fraction below it. This
assessment is report-only: low-Neff variants do not become an active QC-rule
failure and are not removed, but a non-zero count does set the
statistical-quality warning status. The decision is recorded in the terminal
report, canonical log, metric TSV, assessment JSON, run summary, completion
manifest, and HTML report.
See the upstream
[LDSC implementation](https://github.com/bulik/ldsc/blob/master/munge_sumstats.py)
and the general sample-size QC rationale in
[Winkler et al.](https://doi.org/10.1038/nprot.2014.071).

“QC-passed” means that a raw-VCF record passed every active rule under
`modules.qc_summary.rules`. It does not identify another output VCF. The terminal report
shows how many raw records fail each independent condition and explicitly notes
that these counts can overlap. Only the final combined mask determines the
excluded and QC-passed totals.

Persistent reports are written to the configured dataset output tree:

- `<dataset>_chromosomewise_harmonisation_metrics.tsv` contains flattened
  chromosome-wise and dataset-level harmonisation metrics;
- `<dataset>_pre_vcf_column_statistics.tsv` contains per-chromosome statistics
  for the harmonised columns exported for VCF creation;
- `<dataset>_gwas2vcf_column_mapping.json` contains the validated common
  mapping from GWAS-to-VCF field names to exported column positions;
- `<dataset>_<build>_vcf_qc_metrics.tsv` contains raw and final metrics;
- `<dataset>_<build>_qc_summary.csv` combines overall, metric, rule, validation,
  and execution-provenance rows;
- `reports/<dataset>_<build>_qc_report.html` is the detailed human-readable
  view of the same reconciled assessment;
- `<dataset>_<build>_vcf_qc_rule_results.tsv` contains raw-VCF rule criteria,
  actions, and affected-variant counts;
- `<dataset>_<build>_qc_assessment.json` contains the complete structured result.

The temporary extracted TSV is always removed, including when assessment fails.
The raw merged VCF is unchanged and remains the VCF returned by harmonisation.
The bcftools query expressions, scientific rules, report paths, and
temporary-file settings are all declared under `modules.qc_summary`. The
harmonisation resolved-configuration report includes that module because it is
an internal dependency. Screen labels are derived from the resolved VCF mapping,
so a configured tag change is reflected in both analysis and reporting.

`vcf_processing.genome_build_header` configures the explicit VCF metadata line
written during final chromosome concatenation. Its `{build}` placeholder is
resolved separately for each output: the detected input build for source, raw,
and not-lifted records, and the target build for successfully lifted records.
PostGWAS validates that each merged output contains exactly one matching line.
The default is `##genome_build={build}`. It is appended with
[`bcftools annotate --header-line`](https://samtools.github.io/bcftools/bcftools#annotate)
inside the existing uncompressed-BCF merge stream.

### Merged VCF scientific provenance

Every merged VCF also receives the schema-validated metadata declared under
`vcf_processing.provenance`. PostGWAS constructs these values only from the
already resolved sample sheet, study-wide decisions and chromosome resource
preflight. It appends all lines in the same `bcftools annotate` stream used for
the genome-build header and validates every line before accepting the merged
output. No summary-statistics or VCF scan is added.

The default `postgwas_*` headers record:

- the PostGWAS version, dataset ID, source summary-statistics file, full
  resource directory, command-level output directory, dataset-specific output
  directory and run-manifest path;
- the detected source build, the build of that particular VCF, whether
  liftover was applied, and the exact full paths of every resolved FASTA,
  dbSNP VCF, annotation file, chain file and frequency resource used for the
  completed chromosomes;
- the study strand consensus and the exact strand-reference files and
  population-frequency column when that reference was required;
- the EAF source, input column and final EAF/MAF decision, external file
  template and exact chromosome files, allele-swap rule `AF = 1 - input AF`,
  and final `FORMAT/AF` destination;
- the INFO source, input column, external file template and exact chromosome
  files or fixed value, allele-matching rule, interpretation and final
  `FORMAT/SI` destination. An external or fixed INFO value is explicitly
  labelled as a proxy or user-assigned constant, not a study-measured
  imputation-quality statistic;
- the effect, SE and Z input columns, automatic or declared decisions, every
  applied scale/orientation conversion, derivation formula and final
  `FORMAT/ES`, `FORMAT/SE` and `FORMAT/EZ` destinations;
- the p-value input column and detected or declared scale. A supplied raw P is
  retained during harmonisation. A supplied `-log10(P)` is first converted by
  `P = 10^(-input)`. In both cases GWAS-to-VCF subsequently writes
  `FORMAT/LP = -log10(P)`; the two transformations are recorded separately so
  the header never implies that an input `-log10(P)` was passed through
  unchanged;
- trait and case/control inputs, the effective-sample-size rule and
  `FORMAT/NEF` destination, plus the exact population-AF annotation resources
  and configured VCF INFO fields.

In particular, `postgwas_resource_directory` is the resolved
`--resource-directory`, `postgwas_output_directory` is the resolved
command-level `--output-directory`, and
`postgwas_dataset_output_directory` is the exact dataset harmonisation folder
that contains the VCF and its audit outputs.

`postgwas_vcf_created_at` is an ISO-8601 timestamp for creation and structural
validation of that merged VCF. It is intentionally not labelled as overall run
completion, because population-frequency QC and final VCF QC occur afterwards.
The run manifest records the true dataset completion time in `completed_at`
only after those final stages succeed. This distinction preserves an honest
audit timeline while following the field semantics of the
[GWAS-VCF specification](https://github.com/MRCIEU/gwas-vcf-specification).

The supported comparison-frequency panel names and their resource examples are
declared under `comparison_af.available_sources` and
`comparison_af.resource_examples`; the CLI choices and help are generated from
those values rather than maintained independently. This indexed VCF supplies
population INFO annotations and post-merge frequency QC. The packaged default
is ALFA with the EUR population field.

Study AF preflight and the separate external-file classifier use `default_eaf.column`
from the tabular panel selected by `default_eaf.source` and declared by
`resource_layout.default_eaf`. The shipped default is 1000G EUR. The validated
tabular choices are `ALFA`, `wgs_ukb`,
`panukb`, `1000G`, and `fingen`; their exact resource examples are maintained
under `default_eaf.available_sources` and `default_eaf.resource_examples` in
the canonical YAML. Existing run YAMLs that explicitly set `default_eaf.source:
ALFA` retain that selection after updating PostGWAS; change that key to `1000G`
to use the new tabular default. In a full-run YAML the key is under
`modules.harmonisation`. `comparison_af.source` remains the independent VCF
selection, with ALFA as its packaged default.

Every choice uses the same configured chromosome TSV path
contract and must contain the selected population column. This role remains
independent of the comparison VCF annotation panel. Structural columns and the
tab delimiter come from `default_eaf_mapping`. Autosomes are available for all
five configured panels and both supported builds; sex-chromosome coverage is
panel- and build-dependent, so dataset preflight verifies every chromosome
actually required by the study. Direct allele matches use the ALT frequency;
swapped matches use `1-AF`.

For a study-supplied AF column, additional reference comparison is required only
when the whole-study share of usable values at or below 0.5 is **greater than 95%**
(or the configured `eaf.maf_decision_cutoff`). The paired diagnostic set contains
finite, uniquely allele-matched, non-palindromic variants across all chromosomes
at reference sites eligible for biallelic frequency inference. Known extra
alleles exclude a site from these diagnostics; complete site membership in the
supplied reference remains an input assumption. Reference AF is aligned to the
**original study effect allele** before comparing. PostGWAS tests direct EAF
agreement first, without folding either vector to MAF. All configurable
gates below must pass to retain the declared internal EAF source:

| Policy under `modules.harmonisation.policies.eaf` | Default | Acceptance |
|---|---:|---|
| `maf_reference_min_overlap` | 1000 | Matched diagnostic variants ≥ minimum |
| `maf_reference_min_overlap_fraction` | 0.50 | Usable diagnostic pairs / all study rows entering preflight ≥ minimum |
| `maf_reference_min_correlation` | 0.90 | Finite correlation ≥ minimum |
| `study_reference_min_minor_fraction` | 0.90 | Reference aligned effect-allele AF ≤ 0.5 fraction > minimum |
| `maf_reference_max_mean_absolute_difference` | 0.10 | Mean absolute AF difference ≤ maximum |

`maf_reference_correlation_method` remains `pearson` by default; the existing
`spearman` override is also supported. These are dataset QC heuristics, not
validated universal cutoffs or proof that every frequency is correctly labelled.
The usable-overlap denominator is fixed once from all study rows entering AF
preflight after upstream preparation. Diagnostic exclusions (including missing
or invalid frequencies, palindromes, ambiguous matches and non-biallelic sites)
remain in that denominator. Pool usable pairs across all chromosomes; do not
average chromosome correlations or fractions. Exactly 50% passes the default
coverage gate, and the independent 1,000-pair floor still applies. Compare the
fraction directly: 14 usable pairs out of 100 pass a configured minimum of 0.14.
Insufficient count or coverage is inconclusive and cannot select replacement.
The reference-minor numerator includes AF exactly 0.5. Exactly 90% fails its
strict lower bound; 901/1,000 passes. Failure of this gate alone on otherwise
passing direct agreement is also inconclusive and cannot select replacement.

The full run-YAML keys are:

```yaml
modules:
  harmonisation:
    policies:
      eaf:
        maf_reference_min_overlap: 1000
        maf_reference_min_overlap_fraction: 0.50
        maf_reference_min_correlation: 0.90
        maf_reference_max_mean_absolute_difference: 0.10
        study_reference_min_minor_fraction: 0.90
        maf_reference_error_margin: 0.02
```

Proportions use finite values from 0 to 1, not percentages such as `90`.
When every direct acceptance gate passes, that source stays in use even if a
replacement was supplied. Only when coverage is sufficient and direct agreement
itself fails are the same eligible pairs tested for independent MAF support:
folded MAF must pass the overlap, correlation and mean-difference gates, the
reference must pass `reference_minor_fraction_cutoff`, and the MAF error must
beat the EAF error by more than `maf_reference_error_margin`. Only supported
MAF can select the paired [explicit EAF replacement source](eaf_replacement.md).
Its frequencies are a population-reference proxy; the replacement file never
supplies evidence for its own selection or palindrome orientation.

Supported MAF without a replacement, insufficient overlap, constant diagnostic
vectors or otherwise inconclusive evidence stops that dataset before chromosome
analysis, without retrying it. Independent later datasets continue under the
existing batch failure policy; the command ultimately exits nonzero if any
dataset failed. Logs record the measured evidence and failure. Missing and
invalid values remain subject to the existing range/missing policies. These
aggregate gates do not select or remove individual variants, and a mean
difference limit of 0.10 is not a per-variant filter. When replacement is
selected, the original study MAF and source snapshots remain preserved while
the working EAF is obtained from the declared replacement.

When this study-level comparison reads a reference, PostGWAS also retains its
allele matches, including ambiguous and palindromic candidates. Compact matches
are bound to the original source-row IDs in the existing chromosome Parquet work
files. Workers reuse them without rereading the original frequency table or
repeating allele matching or MAF classification. Prepared comparison evidence
does not require an extra CLI option. If the selected replacement is the same
physical reference, population column and structural mapping, workers reuse
its prepared AF instead of reading and staging it again. Other replacements
use the existing external-reference staging and worker memory accounting.
The preparation contract checks reference identity, population, column mappings,
matching policies and unchanged study keys/frequencies. Incomplete or changed
preparation stops instead of silently reusing stale evidence. Preparation is
rebuilt within each invocation; existing validated checkpoint rules still govern
completed-stage reuse. Normal allele transformations and range/missing-value QC
remain in their existing order. Internal matching columns are removed before
scientific outputs. A study that does not trigger the comparison continues to
read its reference normally for chromosome orientation.

The new `study_reference_min_minor_fraction` is a strict **lower bound** for
direct whole-study acceptance. The existing `reference_minor_fraction_cutoff`
remains a separate **upper directional bound** for fallback MAF support; a share
strictly above its default 0.95 makes that classification inconclusive.
If effect alleles are minor, their EAF legitimately equals MAF. The unchanged
`maf_reference_error_margin: 0.02` applies only to fallback MAF classification
after direct EAF agreement fails, and to the separate external-file classifier.
The shared correlation threshold is now 0.90 for direct whole-study agreement,
conditional MAF support and ordinary external-source classification. The existing
0.02 margin remains an aggregate mean-error separation, not a per-variant
tolerance or maximum error; it does not participate in direct acceptance.

A passing direct result means reference-consistent declared EAF, not proof of
the original column's provenance. In particular, a deterministic synthetic
10,000-variant MAF control with 94% minor reference effect alleles produces
direct correlation 0.9523683946, reference-minor fraction 0.94 and direct mean
absolute difference 0.012, and passes all these default gates despite supplying
MAF. These user-selected QC heuristics are not scientifically calibrated to
resolve every EAF/MAF ambiguity. Existing whole-study preparation still requires
full-study data; this change makes no bounded-memory or 100–200-million-row
scalability claim. Allele semantics follow the
[GWAS-SSF contract](https://www.ebi.ac.uk/gwas/docs/summary-statistics-format/);
subsequent swaps follow the
[GWAS Catalog harmonisation rules](https://ebispot.github.io/gwas-sumstats-harmoniser-documentation/Introduction/Harmonising-the-variants/).

The **ordinary sample-sheet external-file classifier** inherits the shared
correlation default of 0.90 and unchanged error margin of 0.02. Its absolute
overlap behavior remains unchanged; the two new whole-study fraction gates do
not apply to it. Its
frequency is screened on its deduplicated raw values before this swap operation. A non-MAF-like result
retains the sample sheet's declared ALT/EAF meaning. A MAF-like result requires
an independent default panel; if `external_eaf_file` and the resource selected
by `modules.harmonisation.default_eaf.source` resolve to the same physical file,
processing stops and asks for a different default source. The frequency column is
accepted as EAF only after at least `eaf.maf_reference_min_overlap`
non-palindromic matches. PostGWAS first folds both frequency columns as
`MAF = min(AF, 1-AF)` and requires the configured Pearson or Spearman
correlation to reach `eaf.maf_reference_min_correlation` and their mean
absolute difference not to exceed
`eaf.maf_reference_max_mean_absolute_difference`. Both requirements must pass:
a high correlation alone can hide a systematic frequency offset, and a missing
correlation from a constant frequency vector is inconclusive. This gate checks
whether the default population panel is compatible enough for the direction
test; folding removes allele direction and therefore cannot prove that the
study field is EAF.

After that suitability gate, the aligned reference must also be informative:
no more than `eaf.reference_minor_fraction_cutoff` of its EAF values may be at
or below 0.5. The column is classified as EAF only when its mean absolute error
as EAF is lower than its error as folded MAF by more than
`eaf.maf_reference_error_margin`; the reverse error separation confirms MAF.
Too few matches, an unsuitable default panel, a reference set dominated by
minor effect alleles, or errors closer than the required margin are
inconclusive. Confirmed MAF and inconclusive results both stop with the measured
evidence; neither is silently used as effect-aligned EAF. Palindromic A/T and
C/G variants do not contribute to this decision. Population-frequency
differences motivating the compatibility gate are documented by the
[1000 Genomes Project Consortium](https://doi.org/10.1038/nature15393); the
resulting EAF must still satisfy the
[GWAS-SSF effect-allele-frequency contract](https://www.ebi.ac.uk/gwas/docs/summary-statistics-format).

All non-VCF reference tables use the same shared, compressed-file-aware reader.
Before a reference reader materializes any rows, PostGWAS checks the delimited
header or Parquet metadata for the configured columns and pushes that
projection into the physical reader. Wide population panels therefore
materialize only chromosome, position, alleles, and the population value needed
by the active analysis.
`external_eaf_mapping`, `external_info_mapping`, and `build_check_mapping`
declare their chromosome, position, first allele, second allele, and delimiter
in the canonical YAML. A delimiter can be explicit or `auto`; automatic
detection is limited to the configured `input.delimiter_candidates` and fails
when no candidate produces a consistent table. The genome-build check paths are
expanded from `resource_layout.build_check` for every build named by
`vcf_processing.target_builds`. Build inference therefore uses configured build
names and reference schemas. It accepts direct, swapped, reverse-complement and
reverse-complement-swapped allele representations but counts each study row at
most once for each candidate build; palindromic alleles, repeated reference
rows, or both reference orientations cannot inflate the evidence. One shared
policy-driven chromosome expression is applied to the study, genome-build
references, strand references, external EAF and INFO tables, concordance
references, and extracted VCF records. It trims labels, optionally removes a
leading `chr`, removes an integral `.0` suffix, removes numeric leading zeroes
when `chromosome.strip_leading_zero` is true, uppercases the result, and then
applies `chromosome.rename_map`. Thus `01`
matches `1` and a configured `23 -> X` mapping is symmetric across the study
and every reference rather than changing the study key alone.
For the separate coverage safeguard, PostGWAS counts unique matched reference
allele records against unique reference allele records only on normalized
chromosome labels present in the study. Input-wide match fractions remain
reported but do not determine the build.
This follows the GWAS Catalog convention of one row per variant while allowing
either other/effect-allele order ([format specification](https://www.ebi.ac.uk/gwas/docs/summary-statistics-format),
[harmonisation methods](https://www.ebi.ac.uk/gwas/docs/methods/summary-statistics)) and accounts for the
documented UCSC/Ensembl chromosome-name difference ([UCSC FAQ](https://genome.ucsc.edu/FAQ/FAQgenes.html)).

Once build inference and study-wide decisions are complete, harmonisation
expands `resource_layout` for only the chromosomes present in that dataset. A
mandatory preflight then validates all resolved files together: VCF/FASTA
indexes and contig compatibility, configured population INFO definitions,
reference-table schemas and data presence, GFF structure, and source-to-target
chain mappings. Shared resources are checked once. Any problem fails the
dataset before chromosome partitions or workers are created, and the complete
grouped problem list is recorded in the dataset log. After that preflight, a
whole-genome user EAF/INFO table is read once into chromosome partitions as
described above; the original source paths, staging row counts, selected
columns, and partition paths are retained in the manifest. No separate
resource-preflight policy is exposed because these checks protect inputs
required unconditionally by the active bcftools/GWAS-to-VCF commands.

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

The audit matches chromosome, minimally represented position, and alleles in
direct, swapped, strand complement, and complement-plus-swapped orientations.
Shared trailing and leading VCF padding is removed only in private comparison
keys; original input and VCF fields remain visible in the reports. This detects
equivalent padded representations without FASTA-backed left-alignment. The
allele-aware match total includes all of those recognized orientations, while
each orientation is also counted explicitly. SNPs, indels, and equal-length
multi-base variants are reported separately. Input-only and VCF-only variants
are descriptive evidence: their number never fails this value-preservation
audit because VCF normalization, multiallelic splitting, and earlier filtering
can change record representation.

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
position counts remain visible. This restrained diagnostic remains useful for
repeat-shifted or otherwise unmatched representations because
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
are written. The default `palindromic_action: compare_resolved` uses each
retained row's `strand_action` from the archived adapter-input audit. This
includes palindromes resolved by internal study EAF when the global consensus
was mixed or unresolved. If that archive is unavailable, validation falls back
to the strong forward/reverse study consensus recorded in the run manifest. No
new palindromic INFO flag is added to delivered VCFs.
The screen describes the recorded `palindromic_comparison_basis`: archived
per-variant actions, study-wide consensus, per-variant actions with consensus
fallback, or mixed recorded evidence. Excluded or unavailable orientation is
labelled explicitly; partially excluded sets retain both counts and the
recorded comparison basis. This wording does not change orientation decisions,
comparison values or the validation verdict.
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
`concordance_validation.staging_batch_rows` bounds the rows held during source-ID
assignment, and `staging_compression` controls the temporary Parquet compression.
Duplicate decisions and a single all-chromosome external EAF file use the same
partition key. A per-chromosome external EAF template is resolved with the VCF's
confirmed build and read only for the active chromosome, avoiding a combined
genome-wide reference table. Each chromosome produces the aligned comparison
frames and counts without calculating metric suites that the combiner would
discard. Compact partition counts are combined globally, and the final overall,
variant-type and position-diagnostic metric suites are calculated once by lazy
scans of the temporary records. Exact medians and Pearson/Spearman correlations
use the combined records rather than averages of chromosome statistics. The
direct Python `compare_input_to_vcf()` interface still returns its full metrics
and verdict. Correlations are unavailable when either compared vector is
constant, because Pearson and Spearman coefficients are undefined in that case
([SciPy Pearson documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.pearsonr.html),
[Spearman documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.spearmanr.html)).
This mathematical precondition adds no tolerance or filtering rule and does
not change counts or validation verdicts. Final reports retain deterministic
source-row ordering, and all staging files are removed on success or failure.

### Consequence annotation chromosome names

`modules.harmonisation.vcf_processing.csq_prefixes_by_build` declares the
literal, case-sensitive prefixes used by native `bcftools csq` in each source
genome build. Every supported build must define `vcf`, `gff` and `fasta`;
`null` means no prefix. The native option order is VCF, GFF, FASTA.

The packaged mappings preserve the prior command's expected layout: VCF and
FASTA names such as `1`, with GFF names such as `chr1`. If all three resources
use matching unprefixed chromosome names, explicitly configure:

```yaml
modules:
  harmonisation:
    vcf_processing:
      csq_prefixes_by_build:
        GRCh37: {vcf: null, gff: null, fasta: null}
        GRCh38: {vcf: null, gff: null, fasta: null}
```

For a `chr`-prefixed GFF with unprefixed VCF/FASTA, set that build's `gff` to
`chr`. This emits `--unify-chr-names -,chr,-`; an all-null mapping omits the
option. Prefixes cannot be empty, contain NUL/commas/whitespace or be the native
`-` sentinel: use YAML `null` for an unchanged label. No CLI flag is added.

Resource preflight scans all GFF feature rows once per unchanged shared file
and compares the configured mapping with the source FASTA index and expected
VCF chromosome labels before workers start. Sequence-region directives alone
do not prove feature coverage. Distinct labels that collapse to one chromosome
are rejected; unrelated alternative contigs need not occur in every resource.
Existing exact FASTA/reference-VCF requirements from normalization and allele
annotation remain authoritative. Prefix matching performs no chromosome alias
inference and does not establish genome-build or sequence equivalence.

The source build, resolved prefixes and expected native lookup labels are
recorded in resource-preflight evidence and annotation command provenance.
Reference FASTA/GFF files and VCF output labels are preserved; no FASTA, GFF or
GTF rewrite is performed. The [`bcftools csq` manual](https://samtools.github.io/bcftools/bcftools.html#csq)
defines the prefix order and preserved output names; the
[native implementation](https://github.com/samtools/bcftools/blob/1.23/csq.c)
requires literal source-prefix removal before adding the destination prefix.
