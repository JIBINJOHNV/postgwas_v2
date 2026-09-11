# Filtering

## Purpose

Filtering materializes a QC-filtered, bgzip-compressed, indexed GWAS-VCF using
one schema-validated filtering configuration. The same resolved configuration
is used by the standalone command and the pipeline.

## What the analysis does

PostGWAS can filter on negative log10 P-value, minor-allele frequency,
imputation quality, study/reference allele-frequency difference, variant type,
palindromic-allele ambiguity, and a build-specific MHC interval. Missing-value
actions are explicit for every nullable field used by an active rule.

`minimum_neglog10_p` is compared with the configured `FORMAT/LP` field. It is
not a raw P-value: genome-wide significance at P ≤ 5×10⁻⁸ corresponds to
approximately `7.30103`.

The genome build is inferred only from the authoritative harmonisation metadata
line `##genome_build=<build>`. Each input must contain exactly one such line,
and its complete value must match a build in `resources.genomes`. Filename text,
`##reference`, contig `assembly`, and liftover metadata are not used as
substitutes. This prevents inherited source-build metadata from overriding the
build explicitly assigned to the merged VCF.

When MHC removal is enabled, filtering selects the matching entry from
`mhc_regions`, requires the header to declare contigs, and requires the
configured chromosome to match a header contig. The equivalent `6`/`chr6`
spelling is accepted and logged. A missing build-specific region or MHC contig
is an error because continuing would silently retain the region.

`--mhc-chrom`, `--mhc-start`, and `--mhc-end` optionally override the selected
build-specific region. Each supplied value applies only to the build inferred
from the current VCF; omitted values continue to come from that build's
`mhc_regions` entry. A coordinate override is rejected unless MHC removal is
enabled, preventing an apparently accepted but inactive region setting. CLI
help displays the GRCh37 and GRCh38 defaults for each option directly from
these schema-validated YAML entries rather than presenting an ambiguous single
coordinate or `unset`.

The packaged one-based, inclusive Genome Reference Consortium MHC intervals are
GRCh37 `6:28,477,797–33,448,354` and GRCh38
`6:28,510,120–33,480,577`. PostGWAS converts the selected start to zero-based,
half-open BED coordinates when it invokes bcftools. These intervals describe
the GRC MHC region, not a broader extended-MHC convention; edit the appropriate
`mhc_regions` entry if an analysis protocol requires a different definition.

## When to use it

Use filtering after harmonisation when you need a reproducible GWAS-VCF subset
for downstream analyses with explicit MAF, INFO, significance, variant-type,
strand-ambiguity, reference-frequency, or MHC rules. Do not use it to repair an
incorrect build, allele orientation, or study-sample mapping.

The default output is a hard-filtered VCF containing only passing variants.
When an audit needs record-level failure evidence, enable
`write_soft_filter_vcf` or pass `--write-soft-filter-vcf`. This writes an
additional VCF that retains every input record and adds one or more semantic
FILTER IDs for the active removal rules that record failed. The hard-filtered
VCF remains the pipeline output used by downstream modules.

## Input requirements

- A non-empty PostGWAS-harmonised `.vcf` or `.vcf.gz` with unique, non-empty
  origin metadata: by default `postgwas_version`, `postgwas_dataset_id`, and
  `postgwas_vcf_status`. Names come from the shared
  `modules.formatting.input_contract.provenance_headers` configuration.
- A VCF header containing exactly one authoritative `##genome_build=<build>`
  declaration, active FORMAT/INFO tags, and,
  when MHC removal is active, contigs.
- The configured `bcftools`, `tabix`, and `bash` executables.
- A filename-safe dataset ID and an output directory.
- The correct reference-population INFO tag when AF concordance is active.

## Direct mode

The standalone entry point is `postgwas sumstat_filter`. The same filtering
service is also available through a configured PostGWAS pipeline.

Filtering resolves bcftools from `resources.executables.bcftools`, whose
packaged value is `bcftools`. There is no filtering-specific `--bcftools`
option; set a nonstandard executable in the run configuration when necessary.

### Hard-filtered VCF

```console
postgwas sumstat_filter \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --dataset-id STUDY \
  --output-directory results/filtering
```

Each filtering switch is enabled by writing the option itself, for example
`--include-indels` or `--write-soft-filter-vcf`. Omitting a switch supplies no
CLI override, so the value continues to come from canonical YAML. These
filtering settings do not expose duplicate `--no-*` aliases or accept a
trailing `true`/`false` value. The packaged filtering defaults for
`include_indels`, `remove_palindromic`, `remove_mhc`, and
`write_soft_filter_vcf` are all `false`, so each behavior is opt-in. A custom
run configuration remains authoritative when its value is not overridden by a
presence flag.

### Explicit filtering policy

```console
postgwas sumstat_filter \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --dataset-id STUDY \
  --output-directory results/filtering \
  --minimum-maf 0.01 \
  --minimum-info 0.7 \
  --minimum-neglog10-p 7.30103 \
  --remove-mhc \
  --mhc-chrom 6 \
  --mhc-start 28477797 \
  --mhc-end 33448354
```

To add the all-record audit VCF without changing the downstream hard-filtered
output, append `--write-soft-filter-vcf`.

Command-line option names are not always identical to the configuration keys
they set: `--minimum-maf` sets `maf_min`, `--minimum-info` sets `info_min`, and
`--maximum-info` sets `info_max`. The common QC-policy options use the same
public names in `sumstat_filter` and the report-only `qc` command, including
`--missing-pvalue-action`, `--missing-af-action`, and
`--missing-info-action`. Use `postgwas sumstat_filter --help` for the current
option names and `postgwas config export --module filtering` for the current
configuration keys.

`info_max` accepts values from 0 through 2, inclusive; its default is `null`,
so no upper INFO filter is applied unless the user selects one. The upper
ceiling supports MaCH Rsq values through 2. Filtering does not infer the score
type: choose the scientifically appropriate `--maximum-info` value for the VCF
being filtered (for example, 1 or 1.05 for a standard INFO policy, or up to 2
for MaCH Rsq).

Filtering and QC summary consume the same validated scientific-policy contract
and have parity regression coverage when the same values are selected. Their
packaged upper-INFO defaults differ: filtering keeps `info_max: null`, whereas
the report-only QC virtual subset uses `info_max: 1.05`. Palindromic and MHC
removal are opt-in in both modules; neither is enabled by the packaged defaults.

## Pipeline mode

Start from an indexed, single-sample PostGWAS-harmonised VCF. Selecting
`sumstat_filter` runs filtering as the requested pipeline analysis:

```console
postgwas pipeline \
  --modules sumstat_filter \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --minimum-maf 0.01 \
  --minimum-info 0.7 \
  --dataset-id STUDY \
  --output-directory results/filtering_pipeline
```

For both the hard-filtered downstream VCF and the additional all-record audit:

```console
postgwas pipeline \
  --modules sumstat_filter \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --minimum-maf 0.01 \
  --minimum-info 0.7 \
  --write-soft-filter-vcf \
  --dataset-id STUDY \
  --output-directory results/filtering_audit_pipeline
```

The pipeline supplies the hard-filtered VCF to subsequent analyses, never the
soft-filtered audit. `--apply-filter` instead requests filtering before another
pipeline target. Some current target combinations, including LD clumping and
fine-mapping, have conflicting MHC CLI options; use a separate filtering run
and supply its validated output to the next pipeline. See
[Pipeline Workflow](../core/pipeline-workflow.md) for the current limitation.

## Parameters

Export the complete canonical configuration before changing less common
policies or field mappings:

```console
postgwas config export --module filtering --style full --output filtering.yaml
```

Important keys include:

- `mhc_regions` (separate GRCh37 and GRCh38 intervals) and the optional
  `mhc_region_override`
- `maf_min`, `info_min`, `info_max`, and `minimum_neglog10_p`
- `missing_af_action`, `missing_info_action`, and `missing_pvalue_action`
- `reference_population_tag` and `frequency_difference_max`
- `include_indels`, `remove_palindromic`, and `remove_mhc`
- `write_soft_filter_vcf` (default `false`)
- `filter_reason_ids`, which defines the stable semantic VCF FILTER IDs
- `vcf_fields`, `output_layout`, `sort_output`, and `report_missing_counts`

VCF tags and output path patterns are validated before execution. In
particular, `reference_population_tag` cannot contain shell syntax, and all
configured outputs must remain below the output directory.

## Processing steps

PostGWAS validates configuration, executables, input metadata, build, contigs,
and active VCF tags before constructing expressions. It then performs a
single-stream reason audit, runs a pipefail-protected bcftools pipeline, creates
the tabix index, counts the result, and reconciles mutually exclusive primary
reasons with the observed number removed.

Immediately after the combined input-validation and record-count stage
completes, the terminal and canonical log show an **Input VCF validation**
block. It records the input basename, authoritative genome build, number of
declared contigs, indexed variant count, and the number of required INFO/FORMAT
fields declared. The successful card also identifies the VCF fields used for
each active scientific measurement and reports configured PostGWAS
field provenance as `AVAILABLE`, `PARTIAL`, or `UNAVAILABLE`. The separate
PostGWAS-origin gate requires the configured version, dataset-ID, and status
headers before variant counting or filtering. Missing, blank, duplicated, or
malformed origin metadata rejects the study VCF and instructs the user to rerun
harmonisation on the original summary statistics, not to add headers manually.
Additional AF, P-value, and score-source metadata remains informational and may
be incomplete. If the field contract fails, the card expands to show every
required field, the rule requiring it, and whether its header declaration is present;
stage 1 then stops before record counting or filtering begins.

Header presence and record-level completeness are distinct. A field reported
as `present` is defined by an appropriate `##INFO` or `##FORMAT` header line;
individual records may still contain `.` for that field. Those record-level
missing-value rules and their removal effects appear in audit order alongside
all other rules in the unified **Filtering results by rule** section.
Related rules are placed beneath field-aware family headings for statistical
significance, allele frequency, imputation quality, study/reference frequency
concordance, variant type, palindromic ambiguity, and genomic region. These
headings only organize the presentation: they never combine or reorder rules.
Every subrule retains its complete name, configured threshold or region, total
failure count, first exclusion count, earlier-rule overlap count, action, FILTER
ID, expression, and machine-readable CSV/TSV representation.

The shared live display reports five measured stages: combined input-contract
validation and input-variant counting, active-rule auditing, hard-filtered VCF
creation and indexing, removal-account reconciliation, and report publication.
A failed stage remains below 100% and is not printed as completed.

After the input card and before rule auditing, the terminal prints a
**Filtering plan** generated from the same ordered rule objects used by the
reason engine. Each row is displayed in the terminal's default colour, numbered,
and marked `EXCLUDE` or `KEEP`. The first failed `EXCLUDE` rule in this displayed
order becomes the variant's mutually exclusive primary removal reason. The final
summary later reports the observed counts under **Filtering results by rule**,
so the plan and its consequences remain visually distinct. The machine-readable
configuration and CSV continue to use the stable internal actions `remove` and
`keep`.

Reason tags are aggregated with exact delimiter-bounded matching in a Polars
streaming query. Thus a configured ID such as `PGWAS_MAF` cannot be confused
with a longer ID that merely contains the same text. The hard non-SNP rule and
its reason tag use the same bcftools `TYPE != 'snp'` expression, including for
mixed multiallelic records. A failed reason audit, unavailable reconciliation,
hard/reason count mismatch, or soft/input record-count mismatch stops the run;
the VCFs and reason report are not published.

The reason-count aggregation runs in a fresh worker with
`POLARS_MAX_THREADS` set before Polars import and verifies that the observed
pool equals resolved `execution.threads`. The requested size, observed Polars
pool, and enforcement status are recorded in the filtering log, summary CSV,
HTML report, and returned result. The same thread value is passed to the final
bcftools BGZF compression stage; bcftools documents that `--threads` primarily
affects compression/decompression rather than the filtering expressions
themselves. See the [Polars thread-pool contract](https://docs.pola.rs/api/python/stable/reference/api/polars.thread_pool_size.html)
and [bcftools scaling guidance](https://samtools.github.io/bcftools/howtos/scaling.html).

Filter expressions and paths are shell-quoted. The VCFs, indexes, MHC BED, and
three report formats are written to temporary paths and published as one
validated artifact set. If publication fails partway through, newly promoted
files are removed and any previous outputs are restored. A failed command
therefore leaves no zero-byte or partially updated successful result; the
canonical log records `STATUS: FAILED` and temporary artifacts are removed.

## Outputs

With the packaged output layout, dataset `STUDY` and build `GRCh37` produce:

- `STUDY_GRCh37_filtered.vcf.gz`
- `STUDY_GRCh37_filtered.vcf.gz.tbi`
- `STUDY_GRCh37_soft_filtered.vcf.gz` and `.tbi` only when
  `write_soft_filter_vcf: true`
- `qc_summary/STUDY_GRCh37_filter_reason_summary.tsv`
- `qc_summary/STUDY_GRCh37_filter_summary.csv`
- `reports/STUDY_GRCh37_filter_report.html`
- `logs/STUDY_GRCh37_filter_gwas_vcf_bcftools.log`
- `STUDY_GRCh37_mhc_exclude.bed` when MHC removal is active

If validation fails before a build can be inferred, the failure is recorded in
`logs/STUDY_filter_preflight.log` instead of fabricating a build-specific name.

## QC and logs

The log records the resolved filtering configuration, resolved executable
paths and versions, active expressions, build/field validation, variant counts,
outputs, runtime, warnings, and final status. An unavailable diagnostic count
is reported as unavailable rather than zero.

The aligned terminal summary, CSV, and self-contained HTML report are rendered
from the same reconciled evidence object. The CSV contains an overall record
and one ordered record for every active rule; its overall row includes input
validation status, input path, declared-contig count, and required-field totals.
The terminal summary ends with the VCF counts before and after filtering and the
retained percentage, placing the overall result after the rule details,
reconciliation, and saved-report paths. The HTML retains its headline outcome
metrics at the top for rapid report navigation.
The HTML begins with a detailed input-validation section containing the header
contract, PostGWAS provenance availability, active measurement fields, and the
complete field-dependency table. It is followed by variant flow, overlap and
first-failure accounting, resolved scientific settings, rule-level counts and
expressions, output links, and execution provenance. Neither report rereads the
VCF or calculates an independent set of scientific metrics.

## Interpretation

The output contains records passing the conjunction of all active include
rules and none of the active exclusion rules. Independent failure counts can
overlap. Primary-reason counts use deterministic first-failure attribution so
they can be reconciled without claiming that a variant failed only one rule.

Every active removal rule reports three related counts:

1. **Variants failing this rule** is the independent number satisfying that
   rule's failure expression. The same variant can be present under multiple
   rules, so these counts must not be summed.
2. **Primary removals assigned here** is the subset first encountered at that
   rule in the recorded audit order. These values are mutually exclusive and
   sum to the number removed from the hard-filtered VCF.
3. **Overlapping earlier removal rules** is the difference between the first
   two values. These variants still fail the current rule and are still removed,
   but their primary removal reason was already assigned to an earlier rule.

The terminal presents that accounting in direct language beneath each rule:

```text
N variants failed this rule and were removed:
    • P failed no earlier EXCLUDE rule
    • O also failed one or more earlier EXCLUDE rules
```

Here, `N = P + O`. The two indented lines distinguish variants first failing at
the current rule from variants that also failed a preceding rule, without
changing the stable audit terminology in the detailed reports.

For each removal rule, the exact identity is:

```text
overlapping earlier removal rules = variants failing this rule - primary removals assigned here
```

The detailed HTML report presents these definitions, the formula, a worked
example selected from the current run, and all three values beside every rule.
The stable CSV and reason-TSV columns `variants_matching_reason` and
`variants_removed_for_this_reason` contain the first and second values,
respectively.

The optional soft VCF retains every input record. Its FILTER column contains
all active removal reasons that apply, so one record can carry several
`PGWAS_*` IDs. A missing-value policy set to `keep` is reported in the reason
summary but is not written as a failure FILTER ID. The reason TSV records the
configured FILTER ID beside each removal rule.

## Common problems

- A missing configured FORMAT or INFO tag stops preflight before filtering.
- A missing, duplicate, or unsupported `##genome_build` declaration stops the run.
- An MHC rule requires the configured chromosome to exist in the VCF header.
- Missing values follow their explicit `missing_*_action`; they are never
  silently treated as passing values.
- An input that already declares a configured `filter_reason_ids` value is
  rejected so an older annotation cannot be misattributed to the current run.

## Limitations

Filtering does not correct upstream allele, build, or metadata errors. It also
does not select one study/sample from a multi-study GWAS-VCF; provide a VCF whose
sample layout already has the intended filtering semantics.
The origin headers are provenance declarations, not cryptographic authentication
or proof of complete scientific QC. They apply to the study input, not external
reference-panel VCFs.

## Scientific references

- [bcftools filtering expressions and target-file behavior](https://samtools.github.io/bcftools/bcftools.html)
- [PLINK 2.0 MaCH Rsq filtering range](https://www.cog-genomics.org/plink/2.0/filter)
- [GWAS-VCF specification](https://github.com/MRCIEU/gwas-vcf-specification)
- [GWAS Catalog summary-statistics format](https://www.ebi.ac.uk/gwas/docs/summary-statistics-format)
- [Genome Reference Consortium GRCh37 MHC interval](https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37)
- [Genome Reference Consortium GRCh38 MHC interval](https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38)
