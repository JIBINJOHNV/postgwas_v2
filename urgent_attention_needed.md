# Urgent Attention Needed

## GRCh37 LD-block resources

The original downloaded LDetect source files are correct: the EUR, AFR, and
ASN populations contain different LD-block boundaries.

The problem occurred during resource preparation. The prepared
`GRCh37_AFR_ldetect.bed.gz` and `GRCh37_EAS_ldetect.bed.gz` files contain the
EUR block coordinates; only their population labels were changed. As a result,
the current annotated VCF does not contain scientifically valid
population-specific AFR and EAS blocks.

Required correction:

- rebuild `GRCh37_AFR_ldetect.bed.gz` from `AFR/fourier_ls-all.bed`;
- rebuild `GRCh37_EAS_ldetect.bed.gz` from `ASN/fourier_ls-all.bed`, mapping
  the upstream ASN name to the PostGWAS EAS population name;
- preserve the original zero-based BED coordinates and generate labels ending
  in `_<START>_<END>`;
- BGZF-compress and index both corrected BED files;
- rerun LD-block annotation and do not use the current annotation output for
  population-specific LD clumping.

Users are responsible for verifying that each supplied population BED has the
intended coordinate provenance and has not accidentally reused another
population's partition.

## Pipeline QC summary can assess the wrong VCF

Status: confirmed pipeline-mode issue; deferred while development focuses on
reviewing and correcting individual modules separately.

This issue does not affect a direct `postgwas qc` run when the user explicitly
provides the original harmonised VCF. It affects a pipeline containing
`qc_summary` when filtering, imputation, or LD-block annotation runs first.

The user provides the original PostGWAS harmonisation VCF and one pipeline
output directory. An earlier pipeline step creates a processed VCF and replaces
the shared `args.vcf` value with that new path. Because `qc_summary` is scheduled
last and resolves its input from the mutated namespace, it assesses the
processed VCF instead of the original user-provided VCF. The separate numbered
QC output directory changes only where reports are written; it does not restore
the correct input.

For example, filtering can reduce 1,000 raw variants to 700. Pipeline QC then
sees 700 variants, can report that all 700 passed and produce a retained
fraction near 100%, while labelling those 700 records the "raw merged VCF". The
real loss of 300 variants is therefore absent from that QC calculation.

Required correction when pipeline integration work resumes:

- preserve the original pipeline-entry VCF as an immutable, named run-context
  artifact before any module executes;
- make the QC runner explicitly assess that original artifact rather than the
  current mutable `args.vcf` value;
- continue allowing downstream analyses to use filtered, imputed, or annotated
  VCF artifacts through explicit context entries;
- record the selected input path, input role, and upstream producer in QC JSON,
  logs, and terminal/tabular summaries;
- choose one canonical owner for raw-merge QC so the harmonisation assessment
  and pipeline assessment do not produce conflicting or duplicate reports; and
- add pipeline regression tests proving that upstream VCF-producing steps
  cannot change the raw VCF assessed by `qc_summary`.

## Require PostGWAS-harmonisation provenance on every downstream VCF

Status: urgent cross-module input-contract requirement; not yet implemented.

PostGWAS harmonisation must add exactly one canonical versioned declaration to
every final VCF intended for downstream analysis:

```text
##postgwas_harmonisation_version=<version>
```

This declaration means that the VCF was created by the PostGWAS harmonisation
module version shown in the value. Harmonisation must continue decomposing
multiallelic records before publishing these VCFs, so downstream modules receive
the biallelic record structure defined by the harmonisation contract.

Required implementation:

- define the declaration name and supported-version policy in canonical,
  schema-validated YAML rather than independently in module code;
- make harmonisation write the declaration to every successful final VCF and
  validate it before publishing the output;
- make every other PostGWAS module and pipeline stage that accepts a VCF require
  exactly one declaration and accept only supported PostGWAS harmonisation
  versions;
- treat a missing, duplicated, malformed, or unsupported declaration as an
  invalid input and fail before analysis or scientific output creation;
- do not infer provenance from the filename, output directory, other VCF
  metadata, or a warning-only fallback;
- centralise the header validation in `postgwas.core` so direct commands and
  pipelines enforce the same contract; and
- add regression tests covering accepted harmonisation outputs and rejection of
  external, unmarked, duplicated, malformed, and unsupported-version VCFs.

## Ambiguous sample size can misclassify a binary study as quantitative

Status: urgent harmonisation-to-formatter scientific risk; deliberately
documented for later correction and not currently being implemented.

The harmonisation sample-sheet generator can map a plain total-sample-size
column such as `N` into the existing control-count field when no explicit
control column is available. With `trait_type: auto`, harmonisation accepts that
lone field as quantitative total N, copies it to effective N, and writes `NCO`
and `NEF` without `NC`. This is correct only when the phenotype is genuinely
quantitative or when the supplied value is already a scientifically valid
effective N for a downstream use that does not require case/control counts.

For a case-control study, total N alone cannot determine effective N, sample
prevalence, or the separate case and control counts required for liability-scale
heritability. For example, 10,000 cases and 90,000 controls have total
`N = 100,000` but balanced-design `Neff = 36,000`. Treating total N as Neff
overstates the effective sample size by 2.78-fold.

The malformed VCF then reaches the formatter with no case counts. The formatter
currently infers quantitative when every `NC` value is missing and binary when
even one `NC` value is present. LDSC formatting therefore writes quantitative
`N`, returns no sample prevalence, and cannot support automatic liability-scale
conversion. Without population prevalence, the heritability module still
produces the ordinary observed-scale result; with population prevalence but no
sample prevalence, it stops before liability-scale analysis. MiXeR and any
other consumer of `NEF` can also receive the overstated value.

Required correction when this issue is scheduled:

- treat `trait_type: auto` with only one control/total-sample-size input as
  scientifically ambiguous and fail before harmonisation creates scientific
  outputs;
- require an explicit quantitative declaration when the lone field is total N,
  and require both positive case and control counts for a binary declaration;
- make sample-sheet generation warn for every automatically selected lone total
  or effective-N column, not only columns recognised specifically as effective
  N, and require the generated draft to be reviewed before execution;
- record the declared trait type and the resolved trait type separately in VCF
  provenance instead of leaving the downstream contract at `auto`;
- add an optional formatter `expected trait type` assertion that validates the
  count-derived result without replacing it, and fails before export when they
  disagree;
- treat partially populated case-count fields as ambiguous rather than allowing
  one non-missing value to classify the complete study silently; and
- add harmonisation, formatter, and heritability regression tests for explicit
  quantitative input, complete binary counts, lone total N under `auto`, one
  stray case count, missing case/control cells, observed-only LDSC, and
  liability-scale LDSC requirements.

LDSC itself does not need a numerical correction: its
[official heritability documentation](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation#conversion-to-liability-scale)
states that observed-scale heritability is the default and that liability
conversion requires both sample and population prevalence. The correction is
to prevent PostGWAS from supplying a scientifically misclassified or
incorrectly sized study.

## Remove `--genome-build` from every module except harmonisation

Status: urgent cross-module CLI and configuration correction; not yet
implemented.

PostGWAS harmonisation is the only module that should ask the user to select a
genome build. Every successfully harmonised VCF already contains exactly one
canonical declaration:

```text
##genome_build=<BUILD>
```

Allowing a downstream module to accept `--genome-build` creates a second source
of truth. A user can select a value that disagrees with the harmonised VCF and
cause the module to choose reference files or genomic intervals from the wrong
assembly. Therefore, remove the public `--genome-build` option from every
PostGWAS module other than harmonisation, including downstream pipeline help and
examples.

Required implementation:

- retain genome-build selection only in the harmonisation CLI and its canonical
  configuration, because harmonisation creates the build-specific VCF;
- make every downstream VCF consumer read and validate the single
  `##genome_build=<BUILD>` declaration with the shared VCF-header contract
  before analysis;
- propagate that validated build through the pipeline run context for later
  modules that consume derived non-VCF artifacts, rather than asking the user
  to declare it again;
- require direct downstream commands that consume only derived non-VCF inputs
  to obtain the build from validated PostGWAS input metadata or a provenance
  manifest, never from a replacement CLI build option;
- remove downstream CLI overrides and module-level build defaults that can
  conflict with the VCF-declared build;
- continue using the validated build internally to select build-specific
  intervals, annotations, LD panels, gene resources, and other references;
- fail before analysis when any selected reference resource does not match the
  VCF-declared build; never silently convert this mismatch into a warning or a
  fallback;
- record the VCF-declared build and every corresponding resource-validation
  decision in the canonical log, checkpoints, and output metadata; and
- add direct-command and pipeline regression tests proving that the build is
  derived from the harmonised VCF, is not user-overridable downstream, and that
  incompatible resources are rejected before scientific output is created.

## Future multi-sample GWAS-VCF handling

Status: future capability; multi-sample VCFs are not part of the current
PostGWAS downstream input contract.

PostGWAS harmonisation currently creates one study/sample per VCF. Until
multi-sample support is deliberately implemented, filtering and every other
downstream VCF module must inspect the VCF sample list and reject zero-sample or
multi-sample inputs before analysis. A dataset ID used for output naming must
not be treated as proof that only that sample will be evaluated.

Required implementation before multi-sample VCFs can be supported:

- add one explicit, schema-validated sample-selection policy shared by direct
  commands and pipelines; never silently aggregate FORMAT values across samples;
- require the requested study/sample to exist in the VCF and fail when sample
  selection is missing or ambiguous;
- subset the selected sample in its own bcftools pipeline stage before applying
  any FORMAT-based include or exclude expression, because sample subsetting and
  filtering in one `bcftools view` invocation can evaluate the filter before the
  sample subset;
- run the scientific filter, filter-reason audit, missing-value accounting,
  variant counts, reconciliation, QC, and reporting on the same selected-sample
  stream;
- publish a clearly identified single-sample output VCF unless a separately
  specified and scientifically validated multi-sample output contract is added;
- preserve and validate PostGWAS harmonisation provenance through sample
  selection; and
- add regression fixtures in which one sample passes and another fails each
  FORMAT rule, plus tests for absent, duplicated, ambiguous, and single-sample
  selections.

## Reinstate explicit VCF dataset-identity validation

Status: urgent provenance safeguard deliberately deferred from the current
release.

The formatter currently requires the configured PostGWAS provenance headers
and records the VCF-declared dataset ID, but it does not compare that value with
the `--dataset-id` used to name new outputs. This permits a user-provided,
single-sample PostGWAS VCF to be processed directly even when its historical
internal dataset/sample name differs from the requested output name. The
formatter selects the dataset declared by the VCF, falling back to its sole
sample when necessary; it does not change variant values.

This temporary behaviour creates a provenance risk: a user can accidentally
provide the wrong study and publish scientifically valid calculations under an
unrelated output identifier. Filename agreement cannot prove study identity.

Required correction in a future release:

- define one schema-validated dataset-identity policy in canonical YAML for all
  downstream VCF consumers;
- distinguish the immutable VCF-declared source identity, selected VCF sample,
  and user-requested output label in the run context, logs, checkpoints, and
  output metadata;
- require exact agreement by default once upstream dataset naming and relabeling
  workflows are available;
- provide an explicit, auditable relabel operation for a verified single-sample
  VCF instead of silently rewriting provenance during formatting;
- reject ambiguous zero-sample and multi-sample inputs before scientific
  extraction; and
- add direct-command and pipeline tests for matching identities, intentional
  relabeling, accidental cross-study input, renamed files, and resumed runs.

## Genome-wide table-processing performance follow-up

Status: future performance work only; reviewed against the current registered
modules and deliberately not being implemented now.

The exact SQLite recommendation applies only to `gcta_cojo`, because it is the
only registered module that currently uses SQLite. On the reviewed synthetic
input of 5 million `.ma` variants plus 5 million BIM variants, its complete
preflight took about 52 seconds at about 50 MB peak memory. The two SQL joins
took about 4.4 seconds; most of the remaining time was Python per-row parsing,
numeric validation, insertion, and maintenance of `PRIMARY KEY` indexes. A
fully vectorised Polars comparison took about 10.8 seconds but about 978 MB,
which is not an acceptable default trade for larger whole-genome inputs.

Required `gcta_cojo` direction when this work is scheduled:

- retain bounded disk-backed validation rather than replacing SQLite wholesale
  with an eager in-memory join;
- replace Python per-row parsing with a projected, vectorised, bounded-batch
  reader while preserving every field-count, finite-value, range, duplicate-ID,
  allele-compatibility, overlap, and line-attribution check;
- benchmark bulk table loading followed by creation of the required unique
  index, and compare it with the existing per-row `PRIMARY KEY` maintenance;
- resolve a scientifically valid completion checkpoint before repeating
  expensive content parsing, while still fingerprinting inputs and validating
  every resumed output under the global resume policy; and
- keep batch sizes, memory limits, temporary storage, and index policy in the
  canonical schema-validated configuration rather than hardcoding them.

The broader rule is generalisable, but SQLite itself is not the shared
solution. Apply bounded, projected and vectorised processing to any operation
that scans millions of variants. Do not replace small gene-, locus-, lead-SNP-
or credible-set loops merely because they use Python.

Highest-priority related module paths:

- `gcta_gene`: the main preflight eagerly reads the complete `.ma` table and
  creates whole-dataset identifier sets, allele dictionaries, VCF-coordinate
  dictionaries, and a complete BIM-ID set. Benchmark a projected streaming or
  disk-backed ID/coordinate/allele join. Its separate GMT pathway conversion
  is already bounded and must not be regressed.
- `magma`: variant preparation eagerly reads the complete p-value and location
  tables into pandas, creates normalised copies, merges them, and optionally
  constructs a whole-dataset ID set for BIM matching. Replace only this
  genome-wide preparation stage with a bounded lazy or disk-backed equivalent;
  preserve lowest-p duplicate selection, input order, coordinates, alleles,
  overlap decisions, and output schemas exactly.
- `imputation` with PRED-LD: chromosome workers return eager frames that are all
  retained in a list before `pl.concat` creates another whole-genome frame.
  Write validated chromosome shards, return paths plus compact metrics, and
  stream-concatenate the final table without retaining every chromosome frame.
- `harmonisation`: whole-genome external EAF/INFO staging, chromosome work, and
  reject concatenation are already bounded. Benchmark spilling the immutable
  initial source snapshot to schema-preserving Parquet so its raw source text
  remains available for reject provenance without pinning those buffers beside
  the transformed dataset throughout the dataset-level stage.
- `formatter`: the current single eager VCF projection is reused across output
  formats, and identical identifier selection is already cached. Change this
  path only if a representative benchmark shows that lower peak memory
  outweighs repeated scans or more complex export orchestration.

No equivalent database redesign is currently indicated for `sumstat_filter`,
`post_imputation_filter`, `annot_ldblock`, `ld_clump`, `finemap`, `heritability`,
`qc_summary`, `mixer`, `pops`, `kpops`, `single_cell`, `magmacovar`, `caldera`,
`flames`, `enrichment`, or `manhattan`. These paths already stream or partition
large inputs, delegate the dominant operation to a validated external tool, or
operate on substantially smaller gene-, locus-, matrix-chunk-, or result-level
data. Reassess them only after module-specific profiling identifies a real
bottleneck.

Scientific and validation requirements for every future optimisation:

- preserve exact retained-variant membership, input/output row accounting,
  duplicate decisions, allele orientation, coordinate matching, null and
  non-finite handling, output ordering, reject reasons, QC metrics, and
  provenance;
- project only required columns, push safe filters before materialisation, and
  use configured batch sizes and bounded temporary storage;
- measure wall time, peak resident memory, and temporary-disk use on both a
  representative 5-million-variant dataset and a full-WGS-scale case;
- add equivalence regression tests comparing all scientific outputs and audit
  counts before and after the optimisation; and
- keep module-specific scientific predicates in their modules while placing
  only genuinely reusable streaming, bulk-load, and bounded-concatenation
  mechanics in `postgwas.core`.
