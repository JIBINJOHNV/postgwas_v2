# LD Clumping

## Purpose

LD clumping identifies approximately independent association signals and
groups them into genomic loci using population-matched linkage disequilibrium
(LD) information.

## What the analysis does

PostGWAS supports one or more configured methods:

- `region` keeps the strongest variant in each population-specific annotated
  LD block.
- `standard` performs two r² clumping passes and merges nearby LD boundaries
  into genomic risk loci.
- `cojo-slct` uses the existing GCTA-COJO stepwise-selection implementation,
  preserves every selected signal and its marginal and joint statistics, and
  groups nearby selected signals into physical loci after model fitting.

The `region` and `standard` methods rank variants by `LP` (-log10 P) rather than
by a recovered P value. For sufficiently large LP, conversion to a Float64 raw
P value underflows to `0.0`. Ranking directly by LP preserves distinctions
between those extreme associations instead of creating artificial ties.

Every selected method is required to succeed. A failed method makes the command
fail, and analysis outputs created or changed during the attempt are renamed
with a `.partial.<run-id>` suffix.

## When to use it

Use LD clumping after harmonisation to identify approximately independent
signals. The `region` method requires LD-block annotations, normally produced
by [LD Annotation](ld-annotation.md). LD clumping is a required pipeline step
for [Fine Mapping](fine-mapping.md). Use reference data whose ancestry and
genome build match the study. Fine-mapping pipelines consume the `standard`
FUMA-style loci; `region` is an alternative block-based summary and its loci
should not be combined with the standard locus definition.

## Input requirements

All methods require an indexed, harmonised, biallelic GWAS-VCF containing the
configured effect, standard-error, allele-frequency, and LP fields. The index
must support `bcftools index -n`, which PostGWAS uses to validate and report the
exact input record count before reference validation. The genome build in the
VCF header must match `modules.ld_clumping.genome_build`.

VCF chromosome, allele, and variant identifiers are read as text before numeric
fields are converted and validated. This preserves contig labels such as `X`,
`Y`, and `MT` even when the schema-inference window contains only autosomes.

The `region` method also requires `%INFO/<POP>_LDblock`. The `standard` method
requires a prepared LD-reference directory containing:

- six files per chromosome that contains a variant passing `lead_pvalue` under
  the default `upper_triangle_dual_index` contract: the forward LD table and
  `.tbi`, reverse LD table and `.tbi`, and allele-aware variant inventory and
  `.tbi`; and
- `ld_reference.yaml`, whose format version, build, population, window,
  minimum r²/MAF, column order, allele-order flag, inventory, and endpoint
  orientation match the run.

The configured missing-chromosome policy determines whether an absent
chromosome file produces a partial-reference warning or a strict error.

Reference contig labels are read from the tabix index, so a reference built
with `chr`-prefixed contigs works without configuration. A chromosome the index
does not contain is reported as an explicit naming mismatch rather than as a
missing index SNP.

Every `region` and `standard` output filename includes the population, so two
populations can share one output directory without overwriting each other.

PLINK 1.9 normally stores each LD pair once. The default
`upper_triangle_dual_index` contract retains that forward table and stores a
reverse tabix sidecar. PostGWAS batches exact index positions against both
files, so a variant is found whether PLINK wrote it as endpoint A or B. The
alternative `symmetric_first_endpoint` contract remains supported and needs
only one symmetrised table.

The `cojo-slct` method requires `--cojo-reference-prefix` naming a matched
PLINK 1 binary prefix with non-empty `.bed`, `.bim`, and `.fam` files. The
reference population is the LD-clumping `--population`, and its genome build is
the LD-clumping `--genome-build`; PostGWAS passes those same resolved values to
the existing GCTA-COJO service. Before formatting the VCF, PostGWAS inspects the
BIM identifier convention so the canonical GCTA `.ma` export uses matching SNP
IDs. GCTA's exact ID, allele, frequency, overlap, sample-count, executable, and
version checks remain unchanged.

## Direct mode

Inspect the current command and accepted options:

```console
postgwas ld_clump --help
```

Export the complete canonical configuration:

```console
postgwas config export \
  --module ld_clumping \
  --style full \
  --output ld_clumping.yaml
```

CLI values override canonical YAML values. Argparse supplies no independent
scientific defaults.

### Region-based pruning

Run annotated-region pruning only:

```console
postgwas ld_clump \
  --vcf study_ldblock.vcf.gz \
  --clumping-methods region \
  --genome-build GRCh37 \
  --population EUR \
  --dataset-id STUDY \
  --output-directory results
```

### Region and standard clumping

Run both packaged methods explicitly:

```console
postgwas ld_clump \
  --vcf study_ldblock.vcf.gz \
  --clumping-methods region standard \
  --ld-folder reference/ld \
  --genome-build GRCh37 \
  --population EUR \
  --dataset-id STUDY \
  --output-directory results
```

### Standard r² clumping

Run standard r² clumping only:

```console
postgwas ld_clump \
  --vcf study.vcf.gz \
  --clumping-methods standard \
  --ld-folder reference/ld \
  --genome-build GRCh37 \
  --population EUR \
  --dataset-id STUDY \
  --output-directory results
```

### GCTA-COJO selection

Run GCTA stepwise selection and group its signals at the configured distance
(250 kb by default):

```console
postgwas ld_clump \
  --vcf study.vcf.gz \
  --clumping-methods cojo-slct \
  --cojo-reference-prefix reference/EUR_plink \
  --genome-build GRCh37 \
  --population EUR \
  --dataset-id STUDY \
  --output-directory results
```

## Pipeline mode

Every recipe starts from an indexed, single-sample PostGWAS-harmonised VCF.
Choose the method explicitly so its required resources and preceding steps are
unambiguous:

| Method | Pipeline prepares | External resource supplied by you |
|---|---|---|
| `region` | LD-block annotation | Population/build-matched BED directory |
| `standard` | No annotation or formatter prerequisite | Prepared indexed pairwise-LD directory with its manifest and inventories |
| `cojo-slct` | GCTA formatter table, with BIM-compatible IDs | PLINK BED/BIM/FAM reference prefix and GCTA executable |

### Region-based pruning

```console
postgwas pipeline \
  --modules ld_clump \
  --clumping-methods region \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --ld-region-dir reference/ld_blocks \
  --ld-block-populations EUR \
  --genome-build GRCh37 \
  --population EUR \
  --dataset-id STUDY \
  --output-directory results/region_pipeline
```

### Standard r² clumping

```console
postgwas pipeline \
  --modules ld_clump \
  --clumping-methods standard \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --ld-folder reference/pairwise_ld \
  --genome-build GRCh37 \
  --population EUR \
  --dataset-id STUDY \
  --output-directory results/standard_pipeline
```

### GCTA-COJO selection

```console
postgwas pipeline \
  --modules ld_clump \
  --clumping-methods cojo-slct \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --cojo-reference-prefix reference/EUR_plink \
  --genome-build GRCh37 \
  --population EUR \
  --dataset-id STUDY \
  --output-directory results/cojo_clumping_pipeline
```

### Region and standard clumping together

```console
postgwas pipeline \
  --modules ld_clump \
  --clumping-methods region standard \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --ld-region-dir reference/ld_blocks \
  --ld-block-populations EUR \
  --ld-folder reference/pairwise_ld \
  --genome-build GRCh37 \
  --population EUR \
  --dataset-id STUDY \
  --output-directory results/combined_clumping_pipeline
```

The BED directory, pairwise-LD directory and PLINK prefix are different resource
contracts, not interchangeable paths. Do not supply the pipeline-generated
annotated VCF or GCTA `.ma` as additional CLI inputs. Keep the outputs labelled
by method: only standard-clumping loci supply the current fine-mapping pipeline.

The current combined `--apply-filter` interface can fail on duplicate MHC
options. Run [Filtering](filtering.md#pipeline-mode) separately first, then use
its validated filtered VCF as `--vcf` here; do not assume help success or a
different resource path repairs that parser conflict.

## Parameters

The packaged configuration follows FUMA-style definitions: independent
significant variants use `lead_pvalue` and `clump_r2`; candidate GWAS variants
are limited by `candidate_pvalue`; lead variants use `lead_r2`; and LD
boundaries within `merge_distance_bp` are merged. Packaged values are P ≤
5×10⁻⁸, candidate P ≤ 0.05, r² 0.6, lead r² 0.1, and 250 kb. All are
schema-validated and configurable.

When `remove_mhc: true`, the configured Genome Reference Consortium coordinates
are applied and the excluded-row count is logged:

- GRCh37: chromosome 6, 28,477,797–33,448,354.
- GRCh38: chromosome 6, 28,510,120–33,480,577.

For `cojo-slct`, PostGWAS derives the MHC exclusion IDs directly from the
matched BIM and passes the resulting configured file to GCTA before stepwise
selection. A selected MHC signal is therefore treated as a failed invariant,
not silently removed after model fitting.

### GCTA-COJO selection and physical loci

GCTA settings remain owned by `modules.gcta_cojo`; LD clumping does not define
a second copy of their defaults. With no overrides, `cojo-slct` therefore uses
GCTA stepwise mode, P ≤ 5×10⁻⁸, a 10,000 kb GCTA model window, 0.9
collinearity cutoff, maximum GWAS/reference A1-frequency difference 0.2,
reference MAF 0.01, and no genomic control. The matching LD-clumping CLI options
override those same canonical GCTA keys.

`--cojo-window-kb` and `--cojo-merge-dist` have different meanings:

- `--cojo-window-kb` is passed to GCTA and changes the joint stepwise model;
- `--cojo-merge-dist` is applied only after GCTA completes and cannot change
  which signals GCTA selected or their `bJ`, `bJ_se`, and `pJ` values.

Genome-wide stepwise selection is run by chromosome when the validated input
supports that execution path. Every worker still receives the complete `.ma`
file for GCTA's phenotypic-variance calculation; `--chr` restricts only the
candidate chromosome. GCTA documents `.jma.cojo` (selected variants and their
joint `bJ`, `bJ_se`, and `pJ`) and `.ldr.cojo` as the outputs of
`--cojo-slct`, while `.cma.cojo` belongs to the separate `--cojo-cond`
analysis. The cross-chromosome LD matrix is block diagonal under GCTA's
documented assumption that variants farther apart than `--cojo-wind` are in
complete linkage equilibrium.

`cojo.parallel_output_contract: selection_only` is the packaged LD-clumping
setting. It merges and validates the chromosome `.jma.cojo` and LD-matrix
outputs, then defines loci without launching a separate full-genome
`--cojo-cond` scan. That scan does not select signals, alter their joint
statistics, or participate in physical-locus grouping. Set the value to
`complete` only when the additional genome-wide conditional `.cma.cojo` table
is an explicitly required deliverable. This setting changes output scope and
runtime, not the COJO selection model. Standalone `postgwas gcta_cojo` retains
the complete output contract by default.

The packaged `--cojo-merge-dist` is 250,000 bp. Signals are sorted by chromosome
and position; a new locus starts when the next selected signal is on another
chromosome or is more than 250,000 bp from the preceding signal. A gap exactly
250,000 bp remains in the same locus. The cited workflow says to use the most
significant index SNP but does not distinguish marginal P from joint `pJ` in
that sentence. PostGWAS makes its interpretation explicit: the packaged locus
label is the member with the smallest marginal GWAS P value. Ties use the
largest absolute marginal Z score, then position and SNP ID.
`cojo.index_pvalue: joint` is an explicit alternative. Both marginal and joint
statistics are retained regardless of the label policy.

For the `region` method, PostGWAS separately audits every variant that passes
`lead_pvalue` but has no selected-population `LDblock` annotation after this
configured MHC exclusion. These variants are outside the annotated blocks
analysed by the region method and therefore cannot become block leads. Missing
annotation does not by itself prove biological LD independence. The terminal
shows the total, while the method log, audit TSV, and HTML retain chromosome,
position, input and canonical IDs, alleles, beta, SE, allele frequency, LP,
reported P value, annotation value, and exclusion reason for every variant.

`missing_index_action: warning_skip` follows FUMA's reference-first logic. An
index absent from the allele-aware inventory, present only with different
alleles, or below `minimum_reference_maf` is excluded with a warning and a row
in the reference-exclusions table. Set the action to `error` for strict mode.
An inventory-verified index with no stored pair or no partner above the active
r² threshold is not missing: it is retained as a self row with r²=1. Following
FUMA, `minimum_reference_maf` is also applied to every LD partner before either
clumping pass; the log records excluded pair-row and distinct-partner counts.

`missing_chromosome_action: warning` applies when a chromosome contains one or
more variants at `lead_pvalue` but its configured LD file or tabix index is
missing or empty. PostGWAS skips only that chromosome, prints and logs a warning,
records the exact excluded significant-variant count and missing paths, and marks
the standard result `partial_reference`; it does not describe that result as
genome-wide. Set the action to `error` to require complete reference coverage.
A chromosome with no variant at `lead_pvalue` needs no LD resource.

The worker pool respects global `execution.threads` and `execution.memory_gb`.
Its memory estimate is controlled by schema-validated
`modules.ld_clumping.compute.minimum_worker_memory_gb` and
`modules.ld_clumping.compute.input_memory_multiplier` values.
The columns shown in the HTML locus, lead-SNP, and independent-significant-SNP
tables come from the schema-validated
`modules.ld_clumping.reporting.result_table_columns` mapping. The packaged
selection keeps those browser tables compact; the linked standard TSV outputs
retain every result and audit column.
Chromosomes with the most significant variants enter the worker queue first.
Each chromosome batches significant positions into one inventory query, one
forward-LD query, and, for upper-triangle data, one reverse-LD query, followed
by one inventory query covering the returned endpoints. Tabix writes those
payloads to temporary files for Polars parsing instead of constructing large
Python strings. Inventory identity checks are vectorised, and the validated
adjacency is reused by both FUMA-style clumping passes. The lead pass tracks
unassigned variants with set membership while preserving the original
descending-LP order. Before the repeated membership joins, PostGWAS restricts
only the GWAS annotation frame to the union of configured windows around
inventory-verified significant indexes. LD-reference partners remain on the
left side of each join and are not removed; GWAS-tagged partners above
`candidate_pvalue` remain identifiable and are excluded after the join exactly
as in the full-chromosome FUMA-style logic.

### Allele handling

The GWAS conversion keeps VCF `REF` as the non-effect allele, VCF `ALT` as the
effect allele, and leaves effect size and effect-allele frequency unchanged.
An allele-order-independent ID is used only for matching GWAS and LD records;
it never changes the scientific allele columns. Duplicate canonical variants
with opposite effect-allele orientation cause a hard error. No automatic allele
flip, beta sign change, or frequency transformation is performed.

### Preparing the standard LD reference

PLINK does not run during `postgwas ld_clump`. It runs only when preparing the
LD reference with `tools/resource_preparation/ld_file_preparation.sh`. The
script calculates allele frequencies once for the requested fileset and runs
the LD calculation for each requested chromosome using explicitly supplied
values. In this workflow, `--keep-allele-order` is always included to preserve
the original A1/A2 encoding of a PLINK 1 binary fileset. Because that option
means A1 is not guaranteed to be the minor allele, the preparation script
validates the reported A1 frequency in `[0, 1]` and records `min(f, 1 - f)` in
the inventory's `minor_allele_frequency` field. It does not reorder A1/A2.
PLINK 2 uses different LD commands and is not a drop-in replacement.

Complete preparation example:

```console
tools/resource_preparation/ld_file_preparation.sh \
  --bfile /path/to/reference-prefix \
  --output-dir /path/to/ld \
  --population EUR \
  --genome-build GRCh37 \
  --chromosomes 1-22 \
  --plink plink \
  --gzip pigz \
  --bgzip bgzip \
  --tabix tabix \
  --window-kb 250 \
  --ld-window-variants 99999 \
  --minimum-r2 0.05 \
  --minimum-maf 0.00001 \
  --threads 8 \
  --chromosome-workers 2 \
  --orientation upper_triangle_dual_index
```

To migrate already indexed PLINK pair files without recomputing LD, add
`--reuse-existing-ld`. The script preserves each forward file, derives the
reverse index, creates an exact allele-aware inventory from the source binary
fileset and PLINK frequencies, and writes the version-2 manifest. Without that
flag it also computes the forward pairs. `symmetric_first_endpoint` can be
selected when a single symmetrised table is preferred. Chromosome work is
bounded by `--chromosome-workers`; the script divides the total `--threads`
budget across workers and reserves per-worker capacity for decompression,
record transformation, sorting, and BGZF compression. Two workers are
appropriate for large references on one SSD because higher concurrency can
increase external-sort I/O contention.

The scientific behavior follows FUMA's `getLD.py`: validate the index against
the reference first, skip an unsupported index, and add a self r² row for a
supported index before applying partner thresholds. PostGWAS strengthens that
logic with exact canonical allele IDs and bidirectional upper-triangle lookup.
See the [FUMA implementation](https://github.com/vufuma/FUMA-webapp/blob/master/scripts/getld/getLD.py)
and [PLINK 1.9 LD documentation](https://www.cog-genomics.org/plink/1.9/ld).

## Processing steps

1. Validate the input GWAS-VCF fields, embedded sample, genome build, and
   indexed record count.
2. Validate the selected reference contracts. For `region`, this is the
   population LD-block INFO declaration in the VCF. For `standard`, this is
   the versioned manifest, build, population, LD window, r²/MAF floors,
   orientation, and allele-order contract. For `cojo-slct`, this is the PLINK
   BED/BIM/FAM set, complete BIM scan and identifier convention, reference
   sample count, and supported GCTA runtime. Standard per-chromosome files are
   checked after the significant chromosomes are known because chromosomes
   without a variant at `lead_pvalue` do not require a reference file.
3. Extract the configured biallelic GWAS fields without changing allele
   orientation or effect statistics.
4. For `region`, audit genome-wide-significant variants without the configured
   population LD-block annotation, select the strongest-LP variant per
   annotated LD block, and write both complete and genome-wide-significant
   pruned results.
5. For `standard`, process chromosomes in a memory-bounded worker pool,
   validate indexes against the inventory, batch both LD endpoint directions,
   and apply the two configured FUMA-style r² clumping stages.
6. For `cojo-slct`, create the canonical GCTA `.ma` file through the existing
   formatter, run the existing GCTA `--cojo-slct` service, merge and validate
   each chromosome's selected signals and joint statistics, and group those
   signals using `cojo.merge_distance_bp`. Under the packaged
   `selection_only` output contract, do not run the scientifically separate
   genome-wide `--cojo-cond` analysis or publish its `.cma.cojo` table.
7. Define standard association boundaries, merge loci within the configured distance,
   renumber loci, and validate method outputs.
8. Render one reconciled terminal summary, machine-readable CSV, and
   self-contained HTML report from the validated in-memory method results, then
   finalise the canonical log. Reporting does not reread the VCF or LD files and
   does not recalculate variants, LD, or loci.

The subordinate `LD clumping analysis progress` display uses the shared
PostGWAS stage-progress system. It reports input GWAS-VCF validation and
method-specific reference validation as separate timed stages, followed by each
selected clumping method and report publication. Their aligned outcomes show
the validated input and reference evidence described above. A failed method
leaves the progress below 100% and identifies the failed stage; report
publication reaches 100% only after both report files have been written
successfully.

Chromosome progress reports excluded significant indexes by their actual
reference-audit reason rather than as an unexplained warning total. Each
reported chromosome shows a GWS input/retained/excluded reconciliation, with
the clumping and exclusion details indented beneath it, for example
`chr1 GWS variants: 1,835 input · 1,825 retained · 10 excluded`, followed by
`Clumping results: 70 independent · 32 lead · 24 loci` and
`Excluded GWS variants: not found 8 · allele mismatch 2`. Missing
chromosome resources, exact-position absence, allele mismatch, and the
reference-MAF floor are counted separately. Any non-exclusion diagnostic
warning is shown separately with a pointer to the detailed log. The final
terminal summary groups all excluded GWS variants beneath their mutually
exclusive LD-reference reason, and the HTML report uses the same terminology.
Both surfaces state the analysis limitation directly: excluded GWS variants
were not used in LD clumping, so additional independent signals or loci may be
missing.

After final locus merging and global locus numbering, PostGWAS compares every
excluded significant index with the inclusive `CHR:START-END` boundaries of
the reported loci. The audit distinguishes
`inside_reported_locus_boundary` from
`outside_all_reported_locus_boundaries` and records an overlapping genomic
locus ID when present. This is coordinate containment only: because the index
failed reference validation, it is never relabelled as an LD member and did not
contribute to clumping, locus boundaries, or locus seeding. Missing-reference
chromosomes necessarily classify outside all reported boundaries because no
locus is created for a skipped chromosome.

## Outputs

Region outputs include the extracted table, all pruned variants, significant
pruned variants, the complete significant-outside-LD-regions audit, and the
region-method log. The audit is written with a header even when its count is
zero. Standard outputs include the
formatted table, method log, genomic-risk-locus summary, hierarchy,
independent-signal boundaries, and a lead-cluster table when loci exist. Exact
paths come from `modules.ld_clumping.output_layout`. The standard method always
writes `standard_reference_exclusions`, including only its header when no index
or chromosome was excluded. For every excluded index, that audit contains its
chromosome, position, canonical and input IDs, exclusion reason, boundary
classification, coordinate-overlapping locus ID, matching locus chromosome,
start and end coordinates, and scientific interpretation.

COJO outputs always include the canonical formatter `.ma`, normalized GCTA
result, GCTA and PostGWAS logs, `cojo_selected_signals`, and `cojo_loci`.
When GCTA selects at least one signal, the raw results and LD matrix are also
retained; `.cma.cojo` is retained only for the `complete` parallel-output
contract, and the exact exclusion list is retained when one is required. The
selected-signal table preserves all normalized GCTA columns and adds locus ID,
index status,
locus bounds, and distance from the preceding selected signal. The locus table
contains the physical bounds, representative index signal, marginal and joint
index P values, member count, and complete selected-signal ID list.

The subordinate formatter and GCTA artifacts use a visible shared hierarchy
with dataset- and population-specific filenames under the requested output
directory:

```text
cojo/
├── formatter/
│   ├── {dataset_id}_{population}_gcta.ma
│   ├── logs/{dataset_id}_{population}_formatter.log
│   ├── reports/{dataset_id}_{population}_formatter_report.html
│   └── run_metadata/{dataset_id}_{population}_formatter_*.yaml
├── raw/
│   ├── {dataset_id}_{population}_slct.jma.cojo
│   ├── {dataset_id}_{population}_slct.ldr.cojo
│   ├── {dataset_id}_{population}_slct.cma.cojo  complete contract only
│   └── {dataset_id}_{population}_slct.log
├── results/{dataset_id}_{population}_slct_gcta_cojo.tsv
├── logs/{dataset_id}_{population}_slct_gcta_cojo.log
├── run_metadata/{dataset_id}_{population}_slct_gcta_cojo_*.yaml
└── chromosomes/               per-chromosome outputs when parallel selection runs
```

The shared directories do not repeat the dataset and population as a directory
level. Instead, every formatter and GCTA artifact filename uses the validated
`{dataset_id}_{population}` subordinate run ID. This keeps complete GCTA
provenance accessible, prevents collisions when an output directory is reused,
and avoids mixing subordinate files with top-level LD-clumping tables and reports.
The paths remain controlled by
`modules.ld_clumping.output_layout.cojo_formatter_directory` and
`modules.ld_clumping.output_layout.cojo_root_directory`; the subordinate naming
patterns are controlled by `modules.ld_clumping.cojo`.

Every successful run also writes:

- `summary_csv`, with overall, per-method, per-chromosome, missing-reference,
  index-exclusion-reason, inside/outside reported-locus-boundary counts,
  individual COJO-warning messages, the configured/observed COJO output scope,
  and output-path records; and
- `html_report`, a self-contained detailed report containing method counts,
  chromosome-level standard-clumping outcomes, reference coverage, runtime MAF
  exclusions, searchable tables containing every genomic risk locus, lead SNP,
  independent significant SNP, COJO physical locus, and COJO-selected signal,
  plus the complete GCTA model settings and warnings, resolved thresholds,
  input/execution provenance, and links to the complete generated artifacts.
  When indexes were excluded, the HTML also contains a searchable variant-level
  table showing which exclusions are outside all reported locus boundaries and
  which overlap a boundary by coordinate only. For an overlapping exclusion,
  the table shows both the matching locus ID and its complete `CHR:START-END`
  boundary; outside-locus rows state that no boundary applies. For the region
  method, it
  contains a separate searchable table with every genome-wide-significant
  variant lacking the configured population LD-block annotation.

The HTML result tables are populated from the validated in-memory frames used
to write the standard and COJO TSVs; the report does not reread or recompute them. Each
table shows every result row, supports text filtering by chromosome, locus,
variant ID, or rsID, and links to its complete TSV for the longer LD-membership
and candidate-list audit columns.

Both filenames are schema-validated, population-scoped patterns under
`modules.ld_clumping.output_layout`. The terminal summary displays their
basenames and uses `logging.terminal_label_width` for one aligned value column.
For COJO, the terminal and HTML report state the exact number of formatted GWAS
variants whose SNP IDs and complete allele pairs match the LD reference, show
that count against the formatted total and percentage, and print every warning
reason and consequence instead of only referring the user to another log.
The CSV, HTML, terminal summary, and canonical stage outcomes all consume the
same result metadata, so reporting cannot create a second scientific result.

The standard result metadata distinguishes `input_significant_variants`,
`significant_variants` actually analysed, `skipped_significant_variants`,
`skipped_chromosomes`, `missing_reference_resources`, and
`reference_coverage_status`. A result with any chromosome or index exclusion
has `status: partial_reference` even when valid loci were produced elsewhere.
If every significant variant is excluded, the analysis fails instead of
publishing an empty result.

## QC and logs

The canonical log and resolved-configuration snapshot record selected methods,
thresholds, genome build, MHC exclusions, reference-manifest properties,
allele policy, per-chromosome counts, stage runtimes and outcomes, commands,
outputs, and completion or failure.

The detailed standard-method log records the reconciled inside/outside counts
and one post-locus-classification line per excluded index. The canonical log
records the same counts and the path to the complete exclusion audit without
embedding a potentially large variant table in checkpoint metadata.

The detailed region-method log records the significant-outside-LD-regions
count, its LP/annotation definition, and one complete line per affected
variant. The canonical log records the count and audit path without embedding
the variant table in run or checkpoint metadata.

Review the input and retained variant counts, genome-wide-significant,
independent, and lead-variant counts, merged-locus count, reference-exclusion
reasons, LD rows before and after the runtime MAF floor, low-MAF partner counts,
and method completion status. A chromosome with no variant passing
`lead_pvalue` does not require an LD file. A missing or empty chromosome LD
file follows `missing_chromosome_action`; an incompatible reference manifest
remains fatal because continuing would mix genome builds, populations, or LD
contracts.

## Interpretation

Lead and independent variants represent association structure under the
selected thresholds and reference panel. A lead variant is not necessarily the
causal variant. Compare the locus summaries with the method-specific log before
using the results for fine-mapping or gene prioritisation.

COJO-selected signals are conditionally independent under GCTA's approximate
joint model and the supplied PLINK panel. They are not FUMA lead SNPs or FUMA
independent significant SNPs. COJO physical loci are a reporting layer over the
unchanged selected set; they are not LD-expanded FUMA boundaries. The current
fine-mapping pipeline continues to consume only the `standard` result.

## Common problems

- The `region` method fails when the selected population's LD-block INFO field
  is absent from the VCF.
- The `standard` method marks results `partial_reference` when a required
  chromosome LD table, reverse sidecar, inventory, or tabix index is absent
  under the default warning policy; strict `error` mode fails instead. A
  missing or incompatible `ld_reference.yaml` is always fatal.
- Low study/reference variant overlap commonly indicates mismatched build,
  ancestry, identifiers, or allele representation.
- Many `allele_mismatch` exclusions usually indicate that the LD panel uses the
  opposite genomic strand; check the panel's build and strand rather than
  flipping study effects automatically.
- Opposite effect-allele orientations for duplicate canonical variants are
  rejected rather than silently flipped.
- `cojo-slct` fails before VCF formatting when the PLINK prefix is missing,
  incomplete, build-incompatible, or uses an unsupported BIM identifier
  convention. A low GWAS/BIM overlap usually indicates an ID or build mismatch.
- A small `.fam` panel is allowed with the existing prominent GCTA warning;
  selected signals should be interpreted cautiously because LD estimates are
  less stable.

## Limitations

Results depend on ancestry matching, genome build, reference coverage,
P/r²/distance thresholds, and the selected method. Region blocks and
FUMA-style standard loci use different definitions and should not be treated
as interchangeable. COJO signals and physical intervals are a third definition
and must also remain labelled separately. LD clumping identifies representative
signals; it does not establish causality.

## Scientific references

- [PLINK 1.9 LD reports and window flags](https://www.cog-genomics.org/plink/1.9/ld)
- [PLINK 1.9 `--keep-allele-order`](https://www.cog-genomics.org/plink/1.9/data)
- [FUMA locus and lead-SNP definitions](https://fuma.ctglab.nl/tutorial)
- [Genome Reference Consortium MHC region](https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38)
- [Berisa and Pickrell LD blocks](https://doi.org/10.1093/bioinformatics/btv546)
- [Official GCTA-COJO options and defaults](https://yanglab.westlake.edu.cn/software/gcta/#cojo)
- [Yang et al. 2012 COJO method](https://doi.org/10.1038/ng.2213)
- [Cross-ancestry meta-analysis interval definition: COJO variants collapsed within 250 kb windows](https://pmc.ncbi.nlm.nih.gov/articles/PMC11383513/)
