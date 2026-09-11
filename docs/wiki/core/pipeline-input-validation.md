# Pipeline input validation

Every downstream pipeline starts from a harmonised PostGWAS GWAS-VCF and uses
the same input-validation session and report format. This layer centralises
evidence, reusable file checks, and presentation. It does not replace the
scientific requirements of the selected methods or silently repair input data.
Standalone harmonisation remains a separate preparation command.

## When checks happen

1. Resolve configuration and CLI overrides, validate the execution plan, and
   report independently missing required arguments.
2. Validate the shared entry VCF and its required index. If this fails,
   dependent module preflights are recorded as blocked.
3. Run each active module's registered resource preflight. Independent module
   failures are collected before refusing analysis; an individual validator may
   stop at its first malformed row or failed dependent requirement.
4. Write the startup audit and display the common validation summary.
5. During execution, validate generated inputs and input-dependent compatibility
   at the existing consuming-stage boundaries. Reuse eligible unchanged resource
   evidence rather than repeating the same file check.
6. Update the audit before announcing successful completion, including observed
   later checks. If a stage raises an exception, save the failure information
   instead. A failed audit write cannot be reported as successful completion.

Required arguments that fail before the validation session begins remain in the
canonical error/screen log; they do not create a file-validation audit. A writable,
safe report destination is required to save the audit.

## Concise screen, complete audit

The terminal uses one separated section per file, with its filename once, the
status of its observed checks, selected file-specific facts and a short check
summary. Standard TBI/CSI/FAI/GZI companions are grouped with their exact parent
path when both have evidence. A companion failure makes the section fail.
Same-named files in different directories are distinguished by their full paths.
No file is reopened to build these sections.

All file families use this shared layout: VCF, PLINK companions, interval tables,
gene and SNP sets, covariate tables, matrices, manifests and opaque resources.
Availability-only resources are labelled as such. Build and population remain
declarations. Contig ranges summarise the set of exact labels; the original
order remains in the audit and no chromosome names are changed.
Warnings, failures and blocked requirements remain visible;
raw metric names and consumers remain in the audit. Pipeline readiness reports
stage readiness rather than repeating the entry VCF's counts and build.
Module fields explicitly identifying a recorded file are combined with that
file's validation evidence. Later findings use its stable file number and show
only new facts; they do not repeat its filename or earlier metadata. Cross-file
findings remain at their consuming stage, referring to these file numbers.
Ambiguous basenames and unattributed messages are not guessed at or discarded.
Native-tool diagnostics, execution commands and analysis-result summaries are
not text-filtered.

## Read the audit by check, not just by status

Audit format version 2 combines each exact `path` into one entry containing
`role`, `roles`, `checks`, `status`, `metrics`, `message`, `messages`, `consumers`
and `observations`. A displayed file also has a `file_number`. Distinct check
observations reference metric keys and message indexes rather than repeat the
path and metric payload. Conflicting values are retained under numbered keys
(for example `rows` and `rows#2`), not overwritten. Availability, structure and
compatibility remain different checks. A cross-file or deferred requirement
may have no single path. The screen counts file sections (with indexes grouped),
not individual observations. A section's status is the most severe
observed status; a successful availability check never overrides a failure.

| Status | Meaning |
|---|---|
| `passed` | The specifically listed checks passed. It does not mean every possible file or scientific check ran. |
| `warning` | A condition is reported under the consumer's existing policy; inspect its meaning before interpreting results. |
| `failed` | The named check or operation failed. An unrelated availability pass does not override it. |
| `blocked` | A prerequisite failed, so the dependent requirement could not be checked. |
| `deferred` | The startup check needs a generated input or a later, input-dependent decision. It is not a pass. |

The audit identifies its `phase` as `pipeline_startup` or `pipeline_execution`.
The latter retains the startup records plus later observed evidence. Its overall
status describes the outcome of that phase, not proof of universal scientific
validity. `stages[].deferred_checks` records the startup boundary: these labels
are not automatically converted to passes, even after successful execution.
Consult the relevant stage's canonical log and scientific output-validation
report for details not yet emitted through the common recorder.

For example, an entry-VCF record can report declared build, sample, contigs and
the count obtained from its index. `bcftools index -n` counts from the CSI/TBI
index; it is not an all-record value validator. Required header declarations
also do not prove that every row contains usable values. [bcftools index manual](https://samtools.github.io/bcftools/bcftools.html#index)

The entry GWAS-VCF must contain exactly one sample column. PostGWAS obtains the
sample list with `bcftools query -l`; zero-sample and multi-sample VCFs fail
before analysis, including a multi-sample VCF that contains the requested
dataset ID. If the only VCF sample ID differs from the run dataset ID, the
single-sample contract remains unambiguous: PostGWAS records a warning, uses
that sole sample for downstream extraction, and asks the user to verify that it
is the intended study. No sample is guessed from a multi-sample input.
[bcftools query manual](https://samtools.github.io/bcftools/bcftools.html#query)

## Current resource coverage

Only resources needed by the selected execution plan and analysis modes are
considered. This table describes current checks, not a promise of exhaustive
content validation for every resource family.

| Resource or consumer | Evidence currently available |
|---|---|
| Entry GWAS-VCF and index | Indexed header, required tag declarations, exactly one sample column, configured build/provenance contract and indexed record count. A lone sample-ID/dataset-ID mismatch is reported as a warning. Row-level values remain the relevant readers' responsibility. |
| PLINK-based consumers | Shared companion-file availability/readability; complete BIM identifier inspection; FAM field counts and sample count, and SNP-major BED dimensions where requested by the consumer. Bundle layout does not establish genotype QC or reference ancestry. |
| LD-block annotation | Exact configured build/population gzip BED4 files are read through EOF; required fields, labels and integer interval bounds are checked on every annotation row. The existing exact contig-name compatibility check is preserved. |
| MAGMA and GCTA-gene resources | Configured gene-coordinate schemas, GMT/native gene-set structure, and fastBAT SNP-set structure/counts. Study/reference overlap and generated gene/set compatibility remain consuming-stage checks. |
| MAGMAcovar and MAGMA cell-type covariates | Existing table checks for fields, unique genes, finite numeric or supported missing values, configured missingness/minimum genes, and nonconstant properties. Compatibility with generated MAGMA results is checked later. |
| PoPS | Feature-matrix dimensions and finite numeric values, gene/feature companion names, configured gene-annotation requirements, and configured control/subset membership. Unchanged exact-contract feature/annotation checks can be reused. |
| K-POPS | Kernel dimensions, finite values and configured symmetry tolerance, with gene-order and annotation requirements. Positive semidefiniteness is not established by the symmetry check. |
| scDRS | Existing H5AD matrix, identifier, annotation and configured expression-filter feasibility checks; native gene-set structure, one-to-one crosswalk parsing and configured cell-covariate requirements. |
| LDSC cell-type analysis | Parsed cell-type manifest contract and the configured reference-file inventory. Availability of an LD-score file does not prove its numeric contents. |
| Standard LD clumping | Reference-manifest schema and declared build, population, orientation and coverage contracts. Candidate-dependent chromosome tables, indexes and variant inventories remain deferred under the configured missing-reference policy. |
| PRED-LD and LDSC heritability | Configured reference/script/merge-alleles/LD-score/SNP-count companion availability. These records explicitly do not establish full content validity. |
| CALDERA and FLAMES | Configured model, script and annotation availability; FLAMES feature-name manifest checks; selected CADD file/index and local VEP-cache availability. Opaque models and cache contents are not certified by these checks. |
| MiXeR | Configured reference/resource/executable availability. Current GO-table checks are limited to the header and first data record, not all rows. |
| Filtering, formatting, plotting and QC | The applicable VCF contracts and existing consumer checks. No unrelated reference panel is required merely because another module supports it. |

BED4 here is a genomic-interval table, not PLINK's binary `.bed`. Interval
coordinates follow the UCSC zero-based, half-open convention; the validator
permits zero-length intervals and does not rewrite coordinates or chromosome
labels. [UCSC BED specification](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)

## BIM identifiers and scientific compatibility

Automatic reference-driven formatter selection scans every BIM row using the
configured field roles, delimiter, rsID pattern, coordinate-ID template and
chromosome-label rules. It reports total, rsID, coordinate-style, unsupported,
missing and duplicate-ID counts. Coordinate-style IDs are compared with the
BIM chromosome, position and allele pair; rsIDs are not independently looked
up in dbSNP by this check.

Mixed or unsupported conventions fail automatic formatter selection instead of
being silently renamed. They are not necessarily invalid PLINK files. Duplicate
counts are reported, while the consuming method retains its full-panel or
matched-subset duplicate policy. Consumers sharing one formatter artifact must
agree on its identifier convention; conflicting reference requirements fail.

BIM's A1/A2 ordering is not guaranteed to mean reference-genome REF/ALT.
PostGWAS recognises either allele-token order when inspecting coordinate-style
IDs, but this is not strand harmonisation or proof of effect orientation.
Genome build and ancestry remain declarations unless a separately reported
method establishes them. [PLINK BIM specification](https://www.cog-genomics.org/plink/1.9/formats#bim)

After formatting, the consumer must still check actual GWAS/reference identifier
overlap, allele compatibility and other method-specific requirements. Knowing
that both files use rsIDs is not enough to prove they contain the same variants.

## Configuration and logging

The report path remains in pipeline YAML; shared presentation defaults are in
the canonical application YAML:

```yaml
pipeline:
  validation:
    report_file: run_metadata/input_validation.yaml
logging:
  file_validation:
    max_screen_files: 20
    max_checks_per_file: 3
    max_list_items: 6
    direct_report_file: 'run_metadata/{command}_input_validation.yaml'
```

`pipeline.validation.report_file` is relative to the configured output root and
must not escape it. The writer refuses symlinks, validated input paths, and
existing files that are not recognised PostGWAS validation audits. It writes
atomically so a partial YAML write is not exposed as the final report.

`logging.file_validation.max_screen_files` defaults to `20`, showing at most 20
ordinary passed file sections and reporting how many additional files remain in
the full audit. Set it to `null` to show every successful file section. The
limit never hides warnings, blocked checks, deferred checks, or failures.
`max_checks_per_file` and `max_list_items` limit check and metric lists on screen,
never warning or failure messages. All raw records remain in the saved audit.
Displayed metric and role labels are configured
in `logging.file_validation.metric_labels` and `role_labels`.
`logging.file_validation.outcome_metric_aliases` maps semantic field labels and
structured stage-detail keys to the existing metric vocabulary; its targets are
schema-validated. A list of targets selects the first metric already recorded
for that file, or the final target if none exists. Thus a count from reading
records is not mislabelled as an index count. These aliases affect presentation,
never which checks run.
The superseded
`pipeline.validation.max_screen_records` key is no longer accepted.
`--hide-screen` suppresses the physical
terminal copy, not the canonical screen transcript or validation audit.
These settings add no scientific thresholds or filtering policies.

## One implementation, one unchanged-input check

The session cache is keyed by exact validator contract and current file
identities. The contract includes the parsing/schema options that affect the
answer. Size, modification/change times, device and inode are checked before
reuse and after a cached operation first runs. This includes file dependencies
inspected by nested cached validators. A pipeline input changed after preflight
is refused. Direct modules may generate a new file version at an earlier path;
that version and changed nested dependencies are revalidated, not reused.
Neither mode carries an old result forward for changed inputs. This cheap in-run check is not
a cryptographic proof and does not replace checkpoint input/output fingerprints.
Another invocation creates a new session and revalidates its inputs.

Modules call the shared functions below. The corresponding local file-check
implementations have been removed; module preflights select paths and pass their
configured requirements. VCF header, sample-list and index-count queries are
shared even when consumers do not explicitly pass a previous result. Explicit
VCF evidence cannot bypass the active session's identity checks.

| Shared implementation in `postgwas.core` | Consumers / file family |
|---|---|
| `vcf` | Pipeline entry/handoffs, filtering, formatting, LD annotation, QC and other VCF consumers: header, exact single-sample validation, indexed count and harmonised-input contracts. |
| `plink`, `variant_identifiers` | Formatter and PLINK-based analyses: companion files, BIM identifiers, FAM structure/counts and BED dimensions. |
| `interval_validation` | LD annotation: gzip BED4 content and exact VCF/BED chromosome-label compatibility. |
| `gene_coordinates` | MAGMA and GCTA-gene coordinate-reference readers, retaining their distinct input schemas. |
| `snp_sets` | GCTA fastBAT set-list structure and counts; bounded, disk-backed duplicate/membership accounting. |
| `gene_sets` | MAGMA GMT/native/column-based membership files and GCTA pathway GMT files: shared structure checks with each consumer's duplicate/empty-set policy. The first MAGMA read also constructs its required records; later membership reads do not repeat static validation. |
| `gene_property_validation` | MAGMAcovar and MAGMA cell-type analysis: native gene-result IDs and covariate tables. |
| `matrix_validation` | PoPS feature chunks/name companions and K-POPS binary kernels. Full matrices are not retained in the cache. |
| `gene_annotation_validation` | PoPS and K-POPS gene/TSS annotation structure. Selected chromosomes, anchor genes and generated-score compatibility remain consumer decisions. |
| `reference_resources` | Reference inventories, contained bundle paths, typed YAML manifests, header/first-record checks and whole-line feature-name lists used by PRED-LD, LDSC, LD clumping, MiXeR, FLAMES and CALDERA. |
| `single_cell_validation`, `ldsc_validation` | scDRS H5AD/cell-covariate/crosswalk/gene-set files and LDSC cell-type manifests, configured LD-score inventory and sumstats headers. |
| `io.tables` | Reusable table readers and required-column checks. Arbitrary full dataframes and mutable analysis outputs are not cached. |

“Once” means **the same file identity and the same validation requirements within
one direct or pipeline invocation**. It does not mean that a file is read only once for all
purposes:

- Cheap file-identity and directory-membership checks still run before reuse.
  They detect replaced files and changed resource bundles without repeating
  matrix scans or table parsing.
- A stricter schema, different policy, or different eligible gene/cell universe
  requires its own check. Existing parsed evidence is reused where it is
  sufficient; a startup pass cannot certify compatibility with data that do not
  exist yet.
- Reading validated data to perform an analysis is not a second validation.
  For example, GCTA still consumes set memberships and matches the actual study
  to its LD reference. These operations must not be removed.
- Newly generated scientific outputs still need validation before downstream
  use. Input evidence cannot certify the output of a transformation.

For example, MAGMAcovar and MAGMA cell-type consumers share the same unchanged
covariate-table check when their policies agree. Once MAGMA generates its gene
results, the covariates must also be checked against that actual gene universe:
missingness and nonconstant-property requirements can differ on that subset.
The shared validator caches compact summaries, not an unbounded gene-by-property
matrix.

Public direct commands reuse the same central read-only checks without enabling
pipeline-only branches or stricter pipeline entry-VCF requirements. Missing or
empty files retain the consuming validator's original policy and error type;
they do not become a generic cache error. An explicitly record-only internal
session still runs every requested check. File sections appear at progress boundaries
(or on exit when there are no internal progress stages). A later new check of
the same file produces only incremental findings; identical summaries are not repeated.
Structured scientific stage outcomes are retained unchanged for downstream
consumers, while their file-specific presentation is combined in the audit.
Detailed direct evidence is
written to `logging.file_validation.direct_report_file`, relative to the output
root. Commands with no recorded checks do not invent a pass or rewrite an earlier
audit; early argument errors remain in the canonical log. Forked workers must
return their evidence explicitly and are not newly certified by this recorder.
Orchestrated checkpoint reuse does not rewrite the audit. Native module
checkpoints finish before the direct audit is published, so they do not own a
partially written version of that audit.
The public direct-command progress boundary includes audit publication, even
when a module also shows its own analysis-stage progress. An audit-write failure
therefore cannot complete the overall command's progress bar.

No scientific threshold is introduced by this presentation change. Native format
invariants
remain enforced, including agreement between the declared binary-kernel byte
width and its reader dtype. A configured scDRS crosswalk is now parsed during
preflight, so malformed crosswalks fail before upstream analysis. Gene-coordinate
integer fields reject non-ASCII digits and Python-only underscore notation rather
than accepting values that are not portable to the external tools' text formats.

Several opaque reference families still receive only availability checks;
header-only checks still do not certify every numeric row. The common report
states that limitation. Existing method-specific scientific checks, generated
output checks and detailed analysis summaries remain; centralising file checks
does not justify weakening or deleting them.

Manhattan R/split-vep probes and FLAMES Python-import probes reuse successful
results for the same executable/script or VCF identities and configured probe
requirements within the session. Runtime-search environment values relevant to
these probes are part of their cache contracts. Do not install or replace runtime
libraries during an active analysis; start a new invocation after environment
changes. FLAMES buffers and replays native probe diagnostics into its module log
when that log becomes available.

## Related pages

- [Pipeline Workflow](pipeline-workflow.md)
- [Input and Output Contracts](input-output-contracts.md)
- [Configuration](configuration.md)
- [Logging and Reproducibility](logging-and-reproducibility.md)
