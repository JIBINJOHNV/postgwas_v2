# How Harmonisation Processes Your Data

PostGWAS takes a study's raw GWAS summary-statistics table and produces
allele-aligned, build-resolved GWAS-VCFs with QC and rejection evidence.

This page is the step-by-step walkthrough of that process. It follows the real
runtime order: the same stage and step names appear in the log files, so you can
read a log alongside this page and know exactly where a run is. Nothing here is
a simplification of a different order.

Read it in three passes if you are new to the command:

1. **What goes in and what comes out** — the boundary of the whole operation.
2. **The complete sequence at a glance** — every numbered step in order, one line
   each.
3. **The stages in detail** — what each step actually does, and why.

## What goes in and what comes out

| You supply | Where it is declared |
|---|---|
| One raw summary-statistics table per dataset | `input_file` in the sample sheet |
| The meaning of each of its columns | the remaining sample-sheet columns |
| A reference-resource tree for the study's genome build | `--resource-directory` or `resources.root` |
| Optional policy overrides | `--run-config` and command-line flags |

| PostGWAS produces | Purpose |
|---|---|
| One merged GWAS-VCF per configured build, bgzipped and indexed | the input to every downstream module |
| Raw GWAS-to-VCF adapter output | audit of the conversion step |
| A rejected-variant file and a reason-by-chromosome matrix | why every discarded row was discarded |
| A field-completeness table | post-parse missingness, configured recovery, and post-recovery final-gate evidence |
| QC reports, an EAF/MAF audit, and a population-frequency comparison | evidence that the file is fit to analyse |
| A run manifest, resolved configuration, command, and logs | reproducibility |
| One row per dataset in `harmonisation_run_summary.csv` | status of every dataset in the run |

A single command drives all of it:

```console
postgwas harmonisation \
  --sample-sheet studies.csv \
  --run-config harmonisation.yaml \
  --resource-directory /absolute/path/to/resources \
  --output-directory results
```

The command first prints a run-level preparation block with the sample sheet,
number of selected datasets, resource directory, output directory, and current
preflight stage. Per-dataset resolved settings follow it; the analysis-start
message is printed only when harmonisation is ready to read that dataset.

While a dataset runs, the subordinate progress display follows the validated
boundaries described below: eight dataset stages, completed chromosomes out of
the observed chromosome total, and five post-merge stages. Parallel chromosome
workers never write to the terminal themselves; the parent advances the shared
counter only after receiving a worker result. The counter reaches 100% only
after every chromosome output and all row-accounting checks required by the
resolved rejection policy validate. Failed chromosomes leave it below 100%
while the configured retry policy is applied. `logging.show_progress` controls
this subordinate detail; the mandatory top-level operation progress and
durable screen transcript remain available.

## The complete sequence at a glance

Seven stages run in this order. Stages 2–6 repeat for every dataset in the
sample sheet; datasets are processed one at a time.

| Stage | Runs | What it achieves |
|---|---|---|
| 1. Run preflight | once per command | Reject an unusable command before reading a large file. |
| 2. Dataset preparation | once per dataset | Clean the study and make the decisions every chromosome must share. |
| 3. Chromosome processing | once per chromosome | Harmonise individual variants and build chromosome VCFs. |
| 4. Failure handling and row reconciliation | once per dataset | Retry, account for every input row, and write rejection provenance. |
| 5. Merge and final QC | once per dataset | Concatenate, annotate, assess, and finish the dataset. |
| 6. Optional input-to-VCF validation | only with `--validate` | Independently re-compare the original input with the merged VCF. |
| 7. Run summary | once per command | Record the outcome of every dataset. |

The 29 numbered steps inside stages 2, 3 and 5 are the ones named in the logs:

| Stage | Log label | What happens |
|---|---|---|
| 2 | Dataset step 1 | Check configuration, column mappings, and policy compatibility. |
| 2 | Dataset step 2 | Check the header and count data rows. |
| 2 | Dataset step 3 | Read and clean the study; reject unusable rows. |
| 2 | Dataset step 4 | Check parsed values and plan missing sample sizes. |
| — | *chromosome-X reconstruction guard* | Stop before further dataset work when X would require assumption-based Z-only BETA/SE reconstruction. |
| 2 | Dataset step 5 | Determine the genome build from the data. |
| 2 | Dataset step 6 | Determine strand policy and consensus. |
| 2 | Dataset step 7 | Determine effect, SE, P-value, and frequency types. |
| — | *exact resource preflight* | Verify every per-chromosome resource for the inferred build. |
| 2 | Dataset step 8 | Split retained study rows and any shared user EAF/INFO table into per-chromosome work files. |
| 3 | Chromosome step 1 | Load the chromosome partition. |
| 3 | Chromosome step 2 | Re-resolve this chromosome's reference files. |
| 3 | Chromosome step 3 | Put effects on the BETA scale. |
| 3 | Chromosome step 4 | Orient alleles against the reference and settle EAF. |
| 3 | Chromosome step 5 | Calculate sample size or effective sample size. |
| 3 | Chromosome step 6 | Recover BETA or SE from Z where possible. |
| 3 | Chromosome step 7 | Harmonise P values. |
| 3 | Chromosome step 8 | Derive a still-missing SE from BETA and P. |
| 3 | Chromosome step 9 | Produce the final Z score. |
| 3 | Chromosome step 10 | Check that BETA, SE, Z, and P agree. |
| 3 | Chromosome step 11 | Obtain INFO. |
| 3 | Chromosome step 12 | Complete variant identifiers. |
| 3 | Chromosome step 13 | Check required output fields and reconcile row counts. |
| 3 | Chromosome step 14 | Export the adapter input table. |
| 3 | Chromosome step 15 | Create the first VCF with the GWAS-to-VCF adapter. |
| 3 | Chromosome step 16 | Normalize, annotate, and lift over with bcftools. |
| 5 | Post-merge step 1 | Merge chromosome VCFs. |
| 5 | Post-merge step 2 | Compare population frequencies. |
| 5 | Post-merge step 3 | Combine side files and clean up. |
| 5 | Post-merge step 4 | Assess the raw merged VCF. |
| 5 | Post-merge step 5 | Validate promised outputs and finish the dataset. |

## How your data is represented at each hand-off

Following the data itself is often the fastest way to understand the order. Each
row below is a real hand-off; the representation changes only where shown.

| After | Your data is | Held as |
|---|---|---|
| Stage 1 | still untouched on disk | your original file |
| Dataset step 3 | one cleaned in-memory table with a stable source-row number | Polars table + raw snapshot for rejections |
| Dataset step 7 | the same table, plus study-wide decisions recorded once | table + resolved metadata |
| Dataset step 8 | retained study rows plus projected shared user EAF/INFO fields | per-chromosome partitions |
| Chromosome step 13 | fully harmonised rows, allele-aligned and statistically complete | per-chromosome table |
| Chromosome step 14 | a TSV plus a JSON column mapping the adapter understands | adapter input |
| Chromosome step 15 | a first, unnormalized VCF in the study's own build | per-chromosome VCF |
| Chromosome step 16 | a normalized, dbSNP- and AF-annotated VCF in both builds | per-chromosome VCF ×2 builds |
| Post-merge step 1 | one whole-genome VCF per build | merged, bgzipped, indexed VCF |
| Post-merge step 5 | a validated, QC-assessed deliverable | the file downstream modules consume |

The rejected rows leave this chain at the step that rejected them and are
written separately with their original values. They never rejoin the retained
path.

## Visual workflow

Read this diagram from top to bottom. The three large boxes show which work is
performed once for the complete study, which work repeats for each chromosome,
and which work happens only after chromosome processing.

![PostGWAS harmonisation workflow showing whole-study, chromosome-wise, and finalization stages](harmonisation-processing-order.svg)

The unresolved-variant branch does not rejoin the retained-variant path. The
retry arrow represents failed chromosomes only; completed chromosomes are not
needlessly repeated.

The most important distinction is:

- **Study-wide decisions are made once.** Genome build, strand consensus,
  beta versus odds ratio, SE scale, P-value representation, INFO versus MaCH
  Rsq representation, and initial EAF/MAF classification must be the same for
  every chromosome.
- **Variant orientation is checked row by row.** After the study-wide decisions
  are known, each variant is compared with its chromosome reference and may be
  kept, swapped, reverse-complemented, rejected, or reported as unresolved.

Therefore, “study strand is forward” does **not** mean PostGWAS blindly leaves
every allele unchanged. It still checks every variant against the reference.

## Key terms

| Term | Plain meaning |
|---|---|
| Effect allele | The allele to which BETA, OR, Z, and EAF refer. |
| Other allele | The second allele in the study file. |
| EAF | Frequency of the effect allele. It changes to `1 − EAF` when effect and other alleles are swapped. |
| MAF | Frequency of the less common allele. It is not necessarily tied to the reported effect allele. |
| INFO | Variant-level imputation-quality score. |
| Neff | Effective sample size, especially for an unbalanced case-control study. |
| REF/ALT | The reference and alternate alleles in a genome-build reference or VCF. |

PostGWAS retains autosomes 1–22 and chromosome X. Before filtering, PLINK 23
is normalized to X, PLINK 25/XY/PAR1/PAR2 to X, PLINK 24 to Y, and PLINK 26/M
to MT. Consequently, pseudo-autosomal variants represented by PLINK code 25
remain on the GRCh X contig, while Y, MT, scaffolds and alternative contigs are
removed and counted. Ordinary Y variants are never relabelled as X. This follows
the [PLINK chromosome-code convention](https://www.cog-genomics.org/plink/2.0/filter)
and [PLINK PAR/VCF behavior](https://www.cog-genomics.org/plink/2.0/data).

The tabular `default_eaf` panel supports allele orientation and MAF/EAF
confirmation. The separate indexed `comparison_af` VCF supplies final VCF
population annotations and frequency QC. Neither silently replaces a missing
study EAF.

## Where each part happens

| Level | How often? | Main purpose |
|---|---:|---|
| Run preflight | Once per command | Check the sample sheet, resource root, configuration, and software. |
| Exact resource preflight | Once per study | Check every resource needed by the inferred build and observed chromosomes. |
| Dataset preparation | Once per study | Clean the input and make decisions shared by all chromosomes. |
| Chromosome processing | Once per chromosome | Harmonise individual variants and create chromosome VCFs. |
| Dataset finalization | Once after chromosome work | Merge VCFs, assess QC, and write final status. |
| Concordance validation | Only with `--validate` | Compare the original input with the same-build merged VCF. |

## Stage 1 — Check the run before reading the large file

**In:** your command, sample sheet, and configuration.
**Out:** a validated run context, or a stop before any large file is read.

PostGWAS performs these checks first:

1. Combine packaged defaults, optional run YAML, and explicit CLI options. CLI
   options have highest priority.
2. Read and validate the version-2 sample sheet. If `--dataset-id` is supplied,
   select only that dataset.
3. Resolve the INFO source. Priority is internal INFO, then external INFO, then
   the explicit `--fixed-info` value when neither file-based source exists.
   With no usable INFO source, the run stops.
4. Check the resource root, output directory, required executables, and the
   configured bcftools liftover plugin. Exact chromosome files cannot yet be
   selected because the study's build and observed chromosomes are not known.
5. Save the resolved configuration, command, normalized sample-sheet row, and
   initial logs. Initialize one screen report per selected dataset and the
   run-level summary CSV with one `NOT_RUN` row for every selected dataset.
6. Process sample-sheet datasets one at a time. A failed dataset does not erase
   another dataset that completed successfully. Before each dataset begins,
   its summary row changes to `RUNNING`.

If one of these checks fails, chromosome analysis has not started.

The screen report is an output-routing step, not a second analysis. PostGWAS
copies each emitted stdout message to the active dataset report and optionally
to the terminal. It performs no additional GWAS scan or scientific calculation.

## Stage 2 — Prepare the complete dataset

**In:** the raw summary-statistics table for one dataset.
**Out:** one work file per chromosome, plus study-wide decisions that every
chromosome will reuse.

These eight steps run once on the complete study.

1. **Dataset step 1 — Check configuration.** PostGWAS checks column mappings,
   alternative inputs, EAF/INFO sources, sample-size declarations, resources,
   and policy compatibility. The full study file is not read yet.

2. **Dataset step 2 — Check the header.** It confirms that every configured
   column exists. It then checks that the input stream is readable and not
   obviously truncated, and counts the data rows for later accounting.

3. **Dataset step 3 — Read and clean the study.** Polars reads the configured
   columns. PostGWAS assigns a stable source-row number and preserves the raw
   values for rejection provenance. It then, in order:

   - reports post-parse missingness for every configured study column, then
     rejects rows missing a read-stage mandatory value or every value in a
     mandatory alternative group such as BETA-or-Z;
   - splits a combined coordinate such as `1:100000`, when configured;
   - normalizes chromosome labels and alleles, and accepts a position only when
     it resolves to a finite whole-number coordinate (`100`, `100.0`, and `1e2`
     all mean position 100; a fractional value is never truncated);
   - rejects missing, malformed, non-integral, non-finite, zero, negative, or
     otherwise sub-minimum coordinates as `invalid_position` under the packaged
     1-based minimum;
   - removes chromosomes outside the configured analysis scope under the
     separate `unsupported_chromosome` reason; the packaged scope retains only
     chromosomes 1–22 and X;
   - converts configured scientific columns to numeric types while keeping
     identifiers as strings; configured missing tokens such as `.` become null,
     and an unparseable scalar numeric cell becomes null without forcing the
     complete column to remain text;
   - rejects missing or non-standard alleles after their canonical uppercase
     representation has been established;
   - when a selected internal INFO column contains values separated by
     `info.multi_value_delimiter`, applies the schema-validated unweighted
     `median` or `mean` aggregation—or stops under `fail`—and records
     invalid-token handling before duplicate ranking; cohort-size weighting is
     never guessed because the input has no aligned cohort-level N vector;
   - identifies variants by chromosome, minimally represented position, and
     normalized unordered allele pair; shared trailing and leading VCF padding
     is removed only in private key columns, while original coordinates and
     alleles remain unchanged and no reference-based left-alignment is done;
     removes and separately reports groups containing both effect/other-allele
     orders; and applies the configured deterministic quality ranking,
     including the resulting scalar INFO, only to scientifically consistent
     same-orientation duplicates.

   The reported input counts are explicitly measured after missing-token
   parsing and before numeric normalization. They therefore do not claim that
   a later numeric conversion has already succeeded. The configured
   `output_layout.field_completeness` TSV and the run manifest distinguish
   each field's input column, selected alternative source, recovery step, and
   final-gate action. After chromosome processing, the same TSV is updated
   with final-gate missing counts and explicit assessed-versus-expected
   chromosome coverage. The packaged `reject` action uses sequential
   attribution in `final_check.require` order; `keep` counts are labelled as
   independent and potentially overlapping, while `fail` is labelled as
   fail-fast. EAF, SE, INFO, SNP, and effective N are not falsely presented as
   mandatory at read time merely because downstream output may need or assess
   them.

   The input ledger is fail-closed and disjoint: variants read must equal rows
   removed for read-stage mandatory values, invalid coordinates, unsupported
   chromosomes, unusable alleles, and duplicates, plus rows ready for
   harmonisation. PostGWAS stops if those counts do not reconcile. Coordinate
   and chromosome-scope validation runs once; it is not repeated after numeric
   normalization because the canonical chromosome and position columns are
   preserved through that conversion.

   Chromosome normalization is not study-only. The exact same resolved
   `chromosome.strip_chr_prefix`, `chromosome.strip_leading_zero`, and
   `chromosome.rename_map` policies are reused for every build, strand,
   external EAF, external INFO and concordance reference before keys are
   compared or deduplicated. For example, study `X` matches reference `23`
   under the packaged map, while study `1` matches reference `01`.

4. **Dataset step 4 — Check the parsed values.** PostGWAS verifies that the
   mapped columns contain usable data. It measures missing sample-size values
   across the whole study. If mean or median filling is configured, one
   dataset-wide fill value is calculated here and reused by every chromosome.
   Too much missing sample size stops the dataset.

5. **Dataset step 5 — Determine the genome build.** Study coordinates and
   alleles are compared with the configured GRCh37 and GRCh38 build-check
   files. Automatic inference requires both strong agreement among positions
   present in either reference and coverage of at least
   `build.min_reference_match_fraction` (default 0.10) of the winning build's
   unique reference markers on normalized chromosomes present in the study.
   Input-wide coverage is logged but is not a decision threshold because the
   build-check files are representative marker panels. Thus chromosome-specific
   input is judged against its chromosome-specific marker panel, while an
   unsupported assembly, incompatible chromosome convention, or corrupt
   reference cannot win from a small matching subset. Zero, weak, or ambiguous
   evidence reaches the same actionable Ambiguous-build failure unless an
   explicit compatible build setting resolves it.

6. **Dataset step 6 — Determine strand policy and consensus.** With
   `strand.mode: auto`, informative non-palindromic SNV matches across the
   complete study are counted as forward or reverse-complement. Palindromic and
   ambiguous matches do not vote. At least 1,000 informative variants and a 99%
   dominant direction are required by default. With
   `strand.mode: reference_aligned`, PostGWAS records the explicit declaration
   instead of inferring consensus, but chromosome rows must still match direct
   or swapped reference REF/ALT.

7. **Dataset step 7 — Determine statistic types.** PostGWAS decides or verifies:

   - BETA versus odds ratio;
   - for odds ratios, log-odds SE versus raw OR-scale SE;
   - complete-study BETA/signed-Z direction agreement for rows that need SE recovery;
   - raw P versus negative-log10 P (other logarithmic representations are refused); and
   - EAF-like versus MAF-like frequency distribution.

   For the P-value column, finite numeric values provide the denominator for a
   study-wide negative-value check. At the packaged 0.1% maximum, any non-zero
   fraction up to and including 0.1% produces a prominent warning; more than 0.1%
   stops before chromosome processing.

   Explicit sample-sheet declarations remain visible and are compared with the
   detected evidence. Automatic beta/OR inference produces a warning because
   an all-positive beta column can resemble odds ratios.

Immediately after dataset step 4, PostGWAS applies the **chromosome-X Z-only reconstruction
guard**. If X is present and the study supplies Z but neither BETA nor SE, the
default `effect_from_z.x_chromosome_z_only_action: fail` stops the dataset
before build/strand inference, resource preflight, partitioning, or worker launch. The autosomal
`2*EAF*(1-EAF)` variance cannot identify male `0/1` versus `0/2` coding,
sex-specific sample composition, or PAR status from pooled EAF and Neff. The
error names three corrections: provide BETA, provide SE, or exclude X by
setting both chromosome allow-lists to chromosomes 1 through 22 only. An
expert may explicitly select `allow_autosomal_assumption`; PostGWAS then warns
and records the assumption in
QC and VCF provenance. Exact recovery from supplied BETA or SE is unaffected.

Before step 8, PostGWAS performs the **exact resource preflight**. It resolves
the source and target build resources for every chromosome actually present in
the study, then checks all problems in one pass:

- required files must exist, be regular files, and be non-empty;
- comparison-AF and dbSNP VCF indexes must be present and readable, and their
  exact chromosome labels must match the study with at least one indexed record;
- the comparison-AF VCF header must define every configured INFO frequency tag
  as an alternate-allele value (`Type=Float`, `Number=A`);
- the source FASTA index must contain the study chromosome label, while target
  FASTA contigs must match the target labels declared by the liftover chain;
- default EAF and configured external EAF/INFO table headers must contain their
  mapped chromosome, position, allele, and value columns, without duplicates,
  and the table must contain data rows; and
- GFF and liftover-chain files must contain a readable record of the expected
  format.

Shared files are checked once. If anything is missing or incompatible,
PostGWAS reports the complete grouped list and starts no chromosome worker.
The resolved per-chromosome resource paths and successful check counts are
stored in the dataset log and run manifest.

The index checks follow the
[bcftools annotation/index requirements](https://samtools.github.io/bcftools/bcftools)
and the standard
[FASTA `.fai` index format](https://www.htslib.org/doc/faidx.html).

8. **Dataset step 8 — Split by chromosome.** This happens only after resource
   preflight passes. Only retained study rows are partitioned. PostGWAS selects
   and writes one chromosome at a time, releasing that subset before selecting
   the next, rather than materialising every chromosome partition together.
   Working rows retain their existing order, and each immutable source subset is
   selected by its stable one-to-one input-row ID. Both are schema- and
   row-count-validated, then atomically published as typed Parquet. The resolved
   `input.chromosome_partition_compression` changes internal storage only, not
   scientific values, types, chromosome assignment, or rejection provenance.
   If a user EAF or INFO path is one file without `{chromosome}` and the study
   contains multiple chromosomes, PostGWAS then scans that file once in
   configured-size batches, projects only the mapped variant and value fields,
   and atomically writes a temporary Parquet partition for each observed
   chromosome. A
   `{chromosome}` path template and a single-chromosome run are left unchanged.
   Staging does not deduplicate, transform alleles, or change annotation values;
   the ordinary chromosome step remains responsible for all scientific
   validation and harmonisation. The full in-memory study is released before
   the external scan. PostGWAS then resolves how many chromosomes may run
   together and how many threads each chromosome receives. The total
   `--threads` budget is shared between adaptive workers (default minimum two
   threads, maximum five). Memory estimates use existing row counts, Parquet
   metadata and reference sizes with configured safety margins. Large and small
   chromosomes are admitted together when their summed estimates fit the
   `--memory-gb` budget after parent headroom; a freed slot is refilled immediately.
   An estimate exceeding the worker budget fails before any chromosome starts.
   `execution.scheduling_mode: fixed` retains fixed per-chromosome reservations.
   Estimates are logged and need calibration; they are not measured peak RAM.
   Before each spawn pool starts, the same
   per-chromosome thread limit is exported through `POLARS_MAX_THREADS`, so
   Polars and bcftools are bounded consistently inside every worker.

At this point PostGWAS has not yet created VCFs. It has a validated study,
recorded study-wide decisions, and one work file per chromosome.

## Stage 3 — Process every chromosome

**In:** one chromosome work file plus that chromosome's reference resources.
**Out:** a normalized, annotated, lifted-over VCF for that chromosome in both
configured builds.

Chromosomes may run in parallel, but the following sixteen steps always occur
in this order inside each chromosome.

1. **Chromosome step 1 — Load the partition.** Read the typed Parquet rows with
   their validated schema and connect them to the raw source snapshot used for
   rejected output. No text re-parse or schema inference occurs here.

2. **Chromosome step 2 — Find the reference files.** Reuse this chromosome's
   already preflighted frequency panel, optional external EAF/INFO, dbSNP,
   source and target FASTA, annotation file, and liftover chain. A lightweight
   existence/index recheck protects against a file being removed between
   dataset preflight and worker startup; the worker does not repeat the full
   schema and bcftools validation.

3. **Chromosome step 3 — Put effects on the BETA scale.** A beta remains beta.
   An odds ratio becomes `BETA = ln(OR)`. If SE was determined to be raw
   OR-scale, it is converted before allele swapping. This ordering is important:
   after conversion, swapping negates BETA while leaving log-odds SE unchanged.

4. **Chromosome step 4 — Orient alleles and obtain EAF.** Every variant is
   matched to the supplied chromosome reference using chromosome, position,
   and alleles. A multiallelic position is not accepted by position alone.
   In automatic mode, SNVs test forward, forward-swapped, reverse-complement,
   and reverse-complement-swapped forms. Reference-aligned mode and all indels
   test direct/swapped allele order
   only; bcftools performs their normalization later. An ordinary
   non-palindromic variant follows its unique row-level match. Before a
   study-supplied EAF is allowed to help orient a palindrome, PostGWAS uses
   those already joined, uniquely oriented non-palindromic rows to compare the
   declared EAF with its `1-EAF` alternative. A non-effect-allele-frequency
   diagnosis requires the configured overlap, negative/inverted correlation,
   maximum inverted error, and error-improvement margin together. A conclusive
   diagnosis stops without changing the data; limited, ancestry-mismatched,
   constant, or near-0.5 evidence remains inconclusive and the user's EAF
   declaration remains authoritative. This adds no reference read or join.

   A palindromic A/T or C/G SNP first follows strong whole-study consensus. When consensus is
   mixed, unresolved, or insufficient, an internal EAF tied to the study effect
   allele may resolve the row only when both frequencies are outside the
   configured ambiguity interval and one orientation passes both the configured
   maximum AF difference and competing-error margin. Strong consensus and
   decisive internal-EAF evidence must agree; a conflict is rejected. PostGWAS
   then updates allele-specific statistics and records the result in
   `strand_action`.

   When no allele ordering matches the population-frequency reference,
   `strand.unmatched_action` rejects the row by default. The explicit `retain`
   choice preserves a row in its supplied effect-allele order; under a strong
   reverse study consensus it first reverse-complements an SNV without swapping
   the effect allele. An unmatched palindrome is eligible for this path only
   when the whole-study consensus is strong. The row is marked as
   `reference_unmatched_retained` (or the reverse-complement form), receives no
   aligned population-reference AF, and must still pass the configured EAF and
   final-completeness policies. Retention does not claim a population-reference
   match. GWAS-to-VCF must subsequently validate its REF allele against the
   build-matched genome FASTA, and exact exported-row versus VCF-record
   accounting stops the chromosome if the adapter skips it. An unmatched A/T
   or C/G palindrome without strong study-wide consensus remains governed by
   `strand.ambiguous_action`; FASTA allele matching alone cannot identify which
   reported study allele direction is correct. `fail` stops the chromosome as
   soon as any unmatched row is detected.

   It then validates the selected internal EAF or, only when no internal EAF
   column was selected, allele-matches the configured external EAF. An external
   EAF file is a dataset-level alternative source; it does not row-wise fill
   null values in a selected study EAF column. If the whole-study screen marked
   an internal column as MAF-like, the reference performs the second check here
   using non-palindromic matches.

   An external source is screened separately on its deduplicated raw frequency
   values before a swapped allele match can apply `1-AF`. A result that is not
   MAF-like retains the sample sheet's declared ALT/EAF meaning. If it is
   MAF-like, PostGWAS requires the same EAF-versus-folded-MAF comparison against
   the independently configured `modules.harmonisation.default_eaf` panel. When
   the external source and that default resolve to the same physical file,
   independent confirmation is impossible, so processing stops and asks for a
   different `modules.harmonisation.default_eaf.source`.

   For either internal or external MAF-like evidence, the
   reference must contain enough matched variants. PostGWAS folds both matched
   columns as `MAF = min(AF, 1-AF)` and requires the YAML-selected Pearson or
   Spearman correlation plus the configured maximum mean absolute difference to
   pass. This checks whether the default population panel is suitable for the
   comparison; because folding removes allele direction, it does not itself
   identify EAF. A constant vector has no finite correlation and is
   inconclusive. The aligned reference must then not be dominated by effect
   alleles whose frequency is already at or below 0.5, and the mean error as EAF
   versus folded MAF must differ by the configured minimum margin. Clear EAF
   evidence for an internal column triggers the deferred
   `1-EAF` correction on allele-swapped rows; clear EAF evidence for an external
   column retains its declared mapping and normal direct/swap alignment. A true
   MAF result, too little overlap, an uninformative reference comparison, or
   errors that are too close stop with an actionable message. PostGWAS never
   guesses that an unresolved MAF-like column is EAF and never auto-converts an
   external MAF column. The population-frequency compatibility rationale follows
   the [1000 Genomes Project](https://doi.org/10.1038/nature15393), while the
   accepted output still follows the
   [GWAS-SSF EAF definition](https://www.ebi.ac.uk/gwas/docs/summary-statistics-format).

   A palindromic frequency from an independent external EAF file cannot orient
   the study beta. After automatic consensus/internal study EAF or an explicit
   reference-aligned declaration has selected the study direction, PostGWAS
   uses the external mapping's declared REF and ALT/effect-allele columns.
   Direct matches use the listed ALT frequency and swapped matches use `1-AF`.
   The panel never changes `strand_action`, BETA, or Z; missing upstream
   orientation is rejected.

   After reference orientation and final effect-frequency alignment are
   complete, PostGWAS repeats duplicate validation on chromosome, position,
   REF and ALT. This catches input rows such as A/G and T/C that were distinct
   before strand resolution but now identify the same variant. Exact agreement
   permits the configured deterministic survivor; a very small configured
   relative tolerance prevents transformation roundoff from creating a false
   numeric disagreement. Conflicting aligned statistics remove the complete
   group. The original rows, step-04 reason and
   per-chromosome duplicate report are retained. With the canonical duplicate
   key, different ALT alleles remain separate. Shared VCF padding is ignored in
   the key, but unnormalised indels are not left-aligned or merged through a
   reverse-complement key.

   Finally, the aligned study EAF is compared with the selected population AF
   as QC. A large difference warns for a non-palindromic row by default, but
   rejects a palindromic row by default because its allele letters cannot
   independently verify the selected orientation.

5. **Chromosome step 5 — Calculate sample size.** Case-control studies use
   cases and controls to calculate
   `Neff = 4 / (1/Ncase + 1/Ncontrol)`. Quantitative studies use total analyzed
   N and cannot enter the case-control formula. The whole-study missing-value
   plan is applied consistently.

6. **Chromosome step 6 — Recover BETA or SE from Z when needed.** Supplied BETA
   and SE are kept. Missing values may be reconstructed from Z, EAF, and Neff
   when scientifically possible. For Z-only reconstruction, EAF must be inside
   `(0,1)`, the effective-variance proxy `2*EAF*(1-EAF)*Neff` must meet its
   configured floor, and the full denominator must be finite and positive.
   See the decision table below.

7. **Chromosome step 7 — Harmonise P values.** Apply the study-wide P-value
   representation and leave one canonical raw p-value column named by
   `pvalue.output_column` (packaged default `PVAL`). Count missing, invalid,
   corrected, or exactly zero values. Before conversion it records whether each
   source value is missing, non-finite, or negative. After either representation
   has produced the canonical raw P value, one shared gate rejects genuine
   missing values as `pval_null` and applies `pvalue.out_of_range` to invalid
   values. A negative `-log10(p)` value is masked during conversion and can
   never become `p=1`. The exact source token and natural-log
   probability are retained privately, so positive raw values below Float64
   and large valid negative-log10 values are not collapsed to a configured
   floor. No canonical negative-log10 column is created here: the GWAS-to-VCF
   adapter consumes raw P text and creates FORMAT/LP only while writing the
   final VCF.

8. **Chromosome step 8 — Derive a still-missing SE.** If SE was neither supplied
   nor recovered from Z, derive it from BETA and P using the configured tail
   assumption. The inverse normal calculation uses `ln(P)` with
   `scipy.special.ndtri_exp`, so probabilities below Float64 retain their
   correct Z and SE. The packaged GWAS default is two-sided. A literal reported
   `p = 0` without SE or usable Z fails by default because zero is
   underflow/censoring, not an exact probability.

9. **Chromosome step 9 — Produce the final Z score.** Keep a supplied Z;
   otherwise calculate `Z = BETA / SE` only where BETA and SE are finite and SE
   exceeds `validation.se_division_floor`. A finite negative BETA is valid, and
   a valid zero BETA gives `Z = 0`. Unusable calculation inputs leave Z empty
   for the single validation gate in the next step; this step removes no row.

10. **Chromosome step 10 — Check statistical agreement.** Apply the six
    basic SE, BETA, and Z checks; handle SE-from-clipped-p-value provenance;
    check `Z ≈ BETA / SE`; compare Z with the harmonised P value using the
    two-sided GWAS relationship and configured tolerance; and apply the
    per-chromosome rejection-fraction guard. The shipped settings warn when this
    step alone removes more than 20% of the variants that entered it and fail
    the chromosome when it removes more than 95%. Rejections from other steps
    are not included in either fraction.

11. **Chromosome step 11 — Obtain INFO.** Use internal INFO first, otherwise
    external allele-matched INFO, otherwise the explicitly requested fixed INFO.
    Missing and out-of-range values follow the resolved policy. A fixed value is
    clearly logged as user-assigned rather than measured.

12. **Chromosome step 12 — Complete variant identifiers.** Preserve and
    normalize a supplied study ID. Fill missing IDs with a deterministic
    chromosome-position-allele ID. dbSNP IDs are assigned later by exact CHROM,
    POS, REF, and ALT matching.

13. **Chromosome step 13 — Check required output fields.** Confirm that every
    required column exists and that required row values are present. A completely
    absent required column fails rather than being skipped.

Before export, PostGWAS verifies:

```text
rows read − rows rejected = rows ready for export
```

An imbalance stops the chromosome because a row was lost or counted twice.

14. **Chromosome step 14 — Export adapter input.** Write the harmonised table,
    JSON column mapping, counts, and summaries required by GWAS-to-VCF. The
    final TSV column, `strand_action`, records the reference-orientation action
    for each retained variant. It is intentionally omitted from the positional
    JSON mapping, so the same TSV is passed directly to GWAS-to-VCF while the
    adapter ignores this audit-only column.

15. **Chromosome step 15 — Create the first VCF.** Run the bundled GWAS-to-VCF
    adapter in the inferred input build. Validation failures return a non-zero
    status. A FASTA access failure stops the chromosome and is reported
    separately from a genuine REF-allele mismatch. PostGWAS then requires the
    raw VCF and index to exist and verifies that its indexed record count is
    exactly the number of harmonised rows exported in step 14. A missing,
    header-only, or partially written VCF therefore cannot be reported as a
    successful conversion. Exact commands, exit status, transcript, and record
    accounting are recorded. This handling follows Python's documented
    [non-zero failure status](https://docs.python.org/3/library/sys.html#sys.exit)
    and pysam's documented
    [FASTA coordinate errors](https://pysam.readthedocs.io/en/stable/api.html#pysam.FastaFile.fetch).

16. **Chromosome step 16 — Normalize, annotate, and lift over.** bcftools:

    1. normalizes against the input-build FASTA and splits multiallelic records;
    2. removes exact duplicate records and clears temporary record IDs;
    3. adds dbSNP IDs using exact allele matching and fills remaining IDs with
       the configured coordinate-allele form;
    4. adds configured population AF fields using exact allele matching;
    5. adds consequence annotations;
    6. reads the annotated VCF header and validates the configured study AF,
       ES, and EZ fields plus the population AF fields discovered from the
       annotation columns;
    7. lifts records to the opposite build while explicitly correcting every
       validated ALT-dependent frequency and signed-effect field;
    8. counts swapped and reference-added records, then either excludes them
       under the shipped `exclude` policy or retains them under `keep`; and
    9. normalizes and splits target-build records, globally sorts them, removes
       exact duplicates, then compresses, indexes, and counts the final VCF.

    Non-lifted records are saved. Plugin rejections follow the configured
    warning/failure thresholds. Policy exclusions, normalization changes, and
    duplicates are reported separately, and a required target VCF with zero
    surviving variants always fails.

### A simple allele-swapping example

Assume one input row contains:

| Field | Study value |
|---|---:|
| Effect allele | A |
| Other allele | G |
| OR | 1.50 |
| EAF | 0.70 |

The chromosome reference says `REF=A` and `ALT=G`. The study effect allele is
therefore the reference allele, but PostGWAS writes the association relative to
ALT. This is a `forward_swapped` variant.

The ordered transformations are:

1. Convert OR before swapping: `ln(1.50) = +0.405`.
2. Set effect allele to G and other allele to A.
3. Negate the effect: `BETA = −0.405`.
4. Invert the frequency: `EAF = 1 − 0.70 = 0.30`.
5. Negate a supplied Z score. Keep a log-odds SE unchanged.

The direction changed because the reported allele changed. The biological
association did not change.

### Possible allele-orientation results

| Result | What PostGWAS does |
|---|---|
| Forward | Keep the direction. |
| Forward swapped | Swap alleles; negate BETA/Z; use `1 − EAF`. |
| Reverse complement | Complement alleles; keep the effect direction. |
| Reverse complement swapped | Complement and swap; negate BETA/Z; use `1 − EAF`. |
| Palindromic resolved by consensus | Use strong full-study non-palindromic strand consensus; decisive internal EAF must agree when available. |
| Palindromic resolved by internal EAF | With no strong consensus, require study and reference AF outside the configured ambiguity interval plus one candidate passing the configured maximum-difference and error-margin rules. |
| Palindromic consensus/frequency conflict | Reject or fail because strong consensus and decisive internal EAF select different study effect directions. |
| Palindromic orientation unavailable | Reject or fail because automatic mode has neither strong consensus nor usable internal study EAF. |
| Palindromic frequency discordant | Reject or fail because both study-EAF candidates exceed the configured reference-AF difference. |
| Palindromic ambiguous | Reject or fail because near-half or non-unique internal study-EAF evidence leaves no safe orientation. |
| Reference-aligned palindrome | Accept one direct/swapped REF/ALT match under the explicit declaration; still reject complement-only or unmatched rows. |
| External palindromic EAF aligned | Use the external mapping's declared ALT AF after upstream orientation; do not change study beta, Z, or `strand_action`. |
| Reference ambiguous | Reject or fail because more than one reference orientation remains. |
| Reference unmatched, default | Reject because no population-reference allele pair matches; `fail` is the strict alternative. |
| Reference unmatched, explicitly retained | Preserve the supplied effect-allele order with no panel AF, or complement an SNV under strong reverse consensus; then require all normal EAF, completeness, and genome-FASTA gates. |
| Unmatched palindrome without strong consensus | Reject or fail as orientation unavailable; `retain` never guesses its study effect direction from FASTA alone. |

### How missing effect statistics are completed

| Values available after effect conversion | Action |
|---|---|
| BETA and a complete SE column | Keep both. Z is kept if supplied, otherwise calculated later. |
| BETA and a partially missing SE column | Preserve every populated SE; for null cells require BETA and signed Z to agree in direction, calculate positive `SE = abs(BETA / Z)`, then use BETA and P as the configured fallback. |
| BETA and Z, but no SE | Require BETA and signed Z to agree in direction, then calculate positive `SE = abs(BETA / Z)`. Opposite signs, or zero BETA with non-zero Z, are rejected by default as `beta_z_sign_discordant`. |
| Z and SE, but no BETA | Calculate `BETA = Z × SE`. |
| Z only | Validate EAF, the configured effective-variance floor, and the finite positive denominator; then estimate standardized BETA and SE using Z, EAF, and Neff. |
| BETA and P, but no SE or Z | Derive SE from BETA and P, then calculate Z. |
| Neither BETA nor Z | Fail because an effect cannot be reconstructed. |

The Z-only calculation produces a standardized effect estimate. It depends on
correct EAF and sample-size assumptions. EAF from a population reference rather
than the measured study sample may make the approximation less accurate. The
default proxy floor is `2*EAF*(1-EAF)*Neff >= 1`; it is a configurable
reconstruction safeguard, not an exact minor-allele-count definition or a
general filter applied to supplied GWAS effects.

Null cells inside an otherwise supplied SE column are first recovered exactly
from sign-consistent BETA and a usable signed Z score. Dataset step 07 checks
that relationship once across the complete study before chromosome fan-out.
Only comparable missing-SE rows enter the denominator, so missing or unusable
statistics and populated study SE cells cannot hide systematic discordance.
The packaged `effect_from_z.max_beta_z_sign_mismatch_fraction: 0.01` allows an
exactly 1% discordant fraction to proceed to row rejection; a strictly larger
fraction stops the dataset before chromosome processing and reports the count,
fraction, columns and configured limit on screen and in the canonical log.

At or below that limit, the packaged
`effect_from_z.beta_z_sign_mismatch: reject` removes and records each row for
which no positive SE can satisfy both supplied statistics. The reject reason is
`beta_z_sign_discordant`, and per-chromosome screen summaries report the count.
The explicit `fail` action remains stricter and stops on any discordant row.
This check applies only where SE is missing and never overwrites a populated
study SE. The packaged
`pvalue.derive_partial_missing_se: true` then applies the configured normal-tail
calculation to any null cells still unresolved. Neither step overwrites a
populated study or Z-derived SE. Setting the policy to `false` disables only the
p-based fallback and leaves the remaining null cells for
`validation.se_invalid`. Literal input `p=0` remains governed separately by
`pvalue.zero_missing_se`.

## Stage 4 — Handle chromosome failures and reconcile rows

**In:** the per-chromosome results, complete or failed.
**Out:** a decided dataset status and one combined rejection file whose row count
is proven against the input.

After each chromosome round:

1. Keep completed chromosomes.
2. Clean partial files and retry failed chromosomes up to the configured limit.
3. Stop retrying when the limit is reached or a retry round makes no progress.
4. Fail the dataset by default if a chromosome still fails. If partial output
   was explicitly allowed, label the result `PARTIAL`, never `OK`.
5. Combine input-stage and chromosome-stage rejected rows into one rejection
   file. The parent streams each shard in `rejects.concat_batch_rows` batches,
   preserves source and row order, and accumulates exact row and reason counts
   per source during the write. It validates the final header and
   empty/non-empty state before deleting source shards; it does not materialise
   all rejected rows or re-read the full gzip output.
6. Fail finalization if required rejection provenance cannot be written.
7. Hard-assert the dataset-wide identity
   `parsed rows = recorded rejects + chromosome exports + unprocessed rows from failed chromosomes`.
   The input-stage reject count must also reconcile with the rows passed to
   chromosome partitioning, every validated partition must be represented by
   exactly one final chromosome status, and completed chromosomes must agree
   with both their partition and reject-shard counts. An `OK` dataset permits
   no unprocessed rows. In an explicitly allowed `PARTIAL` dataset, surviving
   rows from a failed chromosome are reported as unprocessed; they are not
   mislabelled as scientific rejects.
8. Write the reason-by-chromosome report from the same validated streaming
   counts used by the assertion.
9. Finalize the EAF/MAF decision using chromosome reference evidence. This final
   decision replaces the initial MAF suspicion in later QC and concordance.

## Stage 5 — Merge chromosome results and run final QC

**In:** the per-chromosome VCFs.
**Out:** the merged GWAS-VCFs that downstream modules consume, with QC evidence
and a final status.

These five steps run once after chromosome processing.

1. **Post-merge step 1 — Merge chromosome VCFs.** Validate every expected VCF
   and index, create a missing index when possible, and concatenate required
   input-build, target-build, and raw adapter groups. Missing inputs, invalid
   files, failed merges, or missing indexes count as failures. Merges are
   retried, and required merge failure stops the dataset by default. Validated
   inputs are ordered by `chromosome.allowed_after_split` before concatenation;
   build grouping comes from the resolved build pair and configured output
   layout, never from build-like text in the dataset ID. In the same streaming
   merge command, PostGWAS adds and verifies the complete `postgwas_*` VCF
   provenance: input and output directories, exact resource paths, source and
   output builds, strand/EAF/INFO sources, effect/SE/Z/P transformations,
   sample-size derivation, VCF fields, and the VCF creation timestamp. This
   metadata reuses already resolved decisions and adds no variant scan.

2. **Post-merge step 2 — Compare population frequencies.** Query the raw,
   unfiltered, concatenated input-build VCF once. For each configured population,
   report missingness, Pearson correlation, mean absolute AF difference, and
   mean absolute difference after `1-AF`. A sufficiently strong negative
   correlation plus a close and materially improved inverted error is reported
   as `frequency_inversion_suspected` and fails the dataset instead of being
   reduced to an ordinary inconclusive population result.
   Report the closest population and warn when it conflicts with the selected
   population or an external EAF/INFO filename token.

3. **Post-merge step 3 — Combine side files and clean up.** Merge adapter
   summaries and mappings for exactly the chromosomes recorded as completed,
   archive adapter inputs, compress configured outputs, and remove only
   intermediates whose required result was successfully created. Each
   chromosome summary is checked against the canonical adapter-audit schema,
   successful export status, and the independently reconciled chromosome row
   count. The dataset summary is then written atomically with one header. It is
   not reopened and rewritten during QC preparation, so a failed final write
   cannot truncate the previous validated audit.

4. **Post-merge step 4 — Assess the raw merged VCF.** One extraction calculates
   raw VCF metrics and evaluates the configured EAF, INFO, reference-AF,
   variant-type, palindromic, and MHC rules. It also calculates metrics for the
   combined **virtual** QC-passed subset. Rule counts can overlap. The raw VCF
   is not replaced by a filtered VCF.

5. **Post-merge step 5 — Finish the dataset.** Validate all required merged
   VCFs, indexes, chromosome status, and merge status. Write detailed QC, the
   final QC takeaways, combined log, elapsed time, output paths, and terminal
   manifest status. For an `OK` run, remove the raw merged GWAS-to-VCF adapter
   VCF and its exact index by default; `--keep_gwas2vcf_intermediate` or
   `vcf.keep_gwas2vcf_intermediate: true` retains and returns it. The two
   annotated build-specific VCFs are never affected. If finalization stops
   before this successful cleanup boundary, the intermediate can remain for
   diagnosis. Only here does the manifest receive its ISO-8601 `completed_at`
   value; the earlier VCF header timestamp is explicitly the VCF creation time.

In the final QC display:

- ✅ means the summarized requirement passed;
- ⚠️ means the section contains something the user should inspect; and
- ❌ means a required operation failed.

For example, the allele-frequency card shows ⚠️ when even a small number of
variants are discordant or missing reference AF. That warning does not by
itself mean the complete dataset failed.

## Stage 6 — Optional input-to-VCF validation

**In:** the original input table and the merged VCF in the same build.
**Out:** a concordance verdict and per-category difference reports.

This stage runs only with `--validate`, after required VCFs and QC are complete.

1. Select the merged VCF in the inferred input build.
2. Stage the original input, duplicate evidence, and VCF values.
3. Confirm that the VCF header build agrees with the inferred input build.
4. Prepare an optional external EAF source. A single all-chromosome file is
   staged once. A path containing `{chromosome}` is resolved with the confirmed
   build and read only for the chromosome currently being compared.
5. Compare one chromosome partition at a time to limit memory.
6. Report SNP, indel, and other-variant matching separately. SNP orientation can
   consider direct, swapped, complement, and complement-swapped alleles. For all
   sequence alleles, the comparison key removes shared trailing and leading VCF
   padding while preserving the original values. This proves equivalence across
   minimal padding differences but does not use a FASTA or claim equivalence for
   repeat-shifted indels that require left-alignment.
7. Compare oriented effect, final EAF, and Z values using configured tolerances.
8. Write summary, mismatch, input-only, VCF-only, duplicate-VCF, and optional
   all-match reports, then add the result to the run manifest.

A validation failure does not erase harmonised files, but the dataset is not
reported as successfully validated.

The same comparison can be run later against an existing merged VCF with
`postgwas --validate`. See the
[Validation Reference](../reference/validation.md).

## Stage 7 — Finalize the multi-dataset run summary

**In:** the outcome of every dataset in the sample sheet.
**Out:** one auditable row per dataset in the run summary CSV.

After each dataset finishes, fails, or is interrupted, PostGWAS updates that
dataset's row in
`<output>/run_metadata/harmonisation_run_summary.csv`. Its columns follow the
analysis order: dataset identity/status, initial input and ready-row SNP versus
indel/other counts, build inference, study-wide decisions, chromosome evidence,
dataset-wide parsed/rejected/exported/unprocessed accounting, population
comparison, final VCF/QC, and output provenance. A successful row records the
resolved and detected genome build, effect type, P-value type, final frequency
type, population result, key variant counts, and the final manifest. Failed or
partial rows retain the evidence that was available before failure, including
the count of rows explicitly unprocessed because a chromosome failed, and
include an actionable failure reason; datasets not reached after an
interruption remain `NOT_RUN`.

For strand orientation, PostGWAS sums the already-recorded chromosome counters
for completed chromosomes only. It also writes expected, completed, and
summarized chromosome counts and checks:

```text
evaluated variants = retained variants + strand-stage removals
retained variants = forward + forward-swapped
                    + reverse-complement + reverse-complement-swapped
```

The CSV is marked `partial`, `mixed`, or `inconsistent` when chromosome
coverage, status, or arithmetic does not support a complete result. No summary
statistics or VCF is rescanned. After all selected datasets have been handled,
the terminal shows the CSV location and the numbers that are OK, partial,
failed, or not finished.

## What can happen to one input row?

| Final outcome | Meaning | Where to look |
|---|---|---|
| Retained unchanged | The row was already usable and correctly oriented. | Chromosome log and final VCF. |
| Retained after correction | OR conversion, a validated allele swap (including its required `1-EAF` update), or statistic derivation was applied and counted. A suspected non-effect-frequency column is never corrected automatically. | Chromosome log and final VCF. |
| Rejected during harmonisation | An explicit policy removed the row. Original input values and the reason are preserved. | Dataset rejected-variant file and reason report. |
| Unprocessed after chromosome failure | No variant-level rule rejected the row, but its chromosome did not complete. This outcome is permitted only for an explicitly allowed partial dataset. | Dataset manifest reconciliation ledger and run summary. |
| Adapter accounting failure | One or more exported rows were not emitted; the chromosome fails instead of silently continuing. | Adapter transcript and chromosome log. |
| Not lifted | The input-build record exists, but no target-build record was produced. | Non-lifted VCF and liftover counts. |
| Dataset failed | Continuing would be incomplete or scientifically unsafe, or required output/provenance could not be validated. | Dataset log and run manifest. |

## What should the user check before downstream analysis?

1. Final status is `OK` and all expected chromosomes completed.
2. Genome build and strand evidence are convincing.
3. Effect type, SE scale, P-value type, and final EAF/MAF decisions are correct.
4. Rejection counts and reasons are acceptable.
5. Input, exported, adapter, and final VCF counts are understood.
6. Required VCFs and indexes are present and liftover retention is acceptable.
7. Reference-AF concordance and closest-population results are sensible.
8. Raw and virtual QC counts are acceptable for the intended downstream method.
9. Optional concordance passed when `--validate` was requested.

## Related pages

- [Harmonisation Overview](overview.md) — what the command is for and how to run
  it.
- [Harmonisation Sample Sheet](sample-sheet.md) — how to declare your study's
  columns.
- [Harmonisation Outputs and QC](outputs-and-qc.md) — where every file described
  above is written.
- [Harmonisation Configuration](../../modules/harmonisation/configuration.md) —
  the policies and thresholds controlling these decisions.
