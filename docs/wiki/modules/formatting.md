# Formatting

## Purpose

Formatting reads one harmonised GWAS-VCF once and creates validated,
tool-specific inputs for downstream PostGWAS modules. It can also write one
additional user-named table directly from CLI column options.

## What the analysis does

The formatter first verifies the configured PostGWAS provenance, genome-build,
and GWAS-VCF field declarations. It then extracts a typed canonical table with bcftools, selects configured
variant IDs independently for each target, infers quantitative versus binary design from case/control count values,
applies each target's required-field and numeric checks, transforms statistics,
and writes outputs atomically. It can produce MAGMA, GCTA gene, SuSiE, FINEMAP,
PRED-LD, LDSC, and MiXeR inputs and an optional custom table in one run.
For LDSC, an optional HapMap3 merge-alleles table can select variants by both
rsID and allele pair before identifier uniqueness is enforced.

## When to use it

Use it when a downstream analysis requires a tabular contract rather than
GWAS-VCF. Pipeline planning inserts it before dependent analyses and may repeat
it after imputation.

## Input requirements

A PostGWAS-harmonised, biallelic GWAS-VCF; a dataset ID used to name the new
outputs; an output directory; `bcftools`; and
at least one built-in format from CLI/YAML or `--custom-output` with `--id`.
Every configured structural and GWAS-VCF FORMAT field must be declared in the
header; individual values may remain missing where a selected target does not
require them.
`--merge-alleles` is optional for a direct LDSC formatter run and required when
the formatter is planned for the `heritability` workflow.

The formatter is not a second harmonisation stage. The VCF must already contain
the configured PostGWAS version, dataset, status, and genome-build metadata; its
chromosome labels and REF/ALT alleles must already satisfy the configured
canonical patterns. Noncanonical values cause an actionable failure and are
never silently relabelled or uppercased. The numeric cast is still required:
`bcftools query` emits a textual table, so POS and statistical/sample-size
fields must be parsed before downstream calculations and validation. This
contract mirrors the canonical metadata and fields written by the PostGWAS
harmonisation VCF merge stage.

For the current release, the formatter records the dataset declared in VCF
provenance but does not require it to equal `--dataset-id`. The requested
dataset ID controls output naming only. The VCF must contain exactly one sample
column. A lone sample ID that differs from the run dataset ID is used with a
warning; zero-sample and multi-sample VCFs fail before extraction.
Dataset-identity validation is deferred in `urgent_attention_needed.md`.

## Command

```console
postgwas formatter --vcf PATH --output-directory PATH \
  [--format FORMAT [FORMAT ...]] [--custom-output FILE] [options]
```

Standalone formatter runs require an explicit `--output-directory`; the
formatter does not silently select a destination. Pipeline runs continue to
use the pipeline's required output directory and create a formatter step
subdirectory within it.

## Minimal example

```console
postgwas formatter \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --dataset-id STUDY \
  --output-directory results \
  --format ldsc \
  --merge-alleles reference/w_hm3.snplist
```

Create an editable formatter YAML containing only shared settings and the LDSC
target-specific sections:

```console
postgwas config export \
  --module formatting \
  --format ldsc \
  --style minimal \
  --output formatting.yaml
```

Multiple target names are accepted after `--format` and are written in the
configured formatter order. Without `--format`, configuration export retains
all formatter target schemas. The selector is for `--module formatting` only;
pipeline export continues to derive required formatter targets from the
pipeline plan.

### Custom CLI output example

The custom table is additional to any built-in formats and requires no user
YAML changes. The output columns follow the order in which their options appear
on the command line.

```console
postgwas formatter \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --dataset-id STUDY \
  --output-directory results \
  --format magma \
  --custom-output STUDY_for_my_tool.tsv \
  --id SNP \
  --chr CHR \
  --pos BP \
  --alt A1 \
  --ref A2 \
  --beta BETA \
  --se SE \
  --p P \
  --eaf EAF \
  --n N \
  --variant-id-type unique
```

## Full example

```console
postgwas formatter \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --dataset-id STUDY \
  --output-directory results \
  --format magma gcta_gene susie finemap pred_ld ldsc mixer \
  --duplicate-id-policy exclude_all \
  --run-config formatting.yaml \
  --resume
```

The formatter resolves `resources.executables.bcftools` (default `bcftools`)
from the run configuration and validates that command on `PATH` before reading
the VCF.

## Parameters

`--format` accepts exactly `magma`, `gcta_gene`, `susie`, `finemap`, `pred_ld`,
`ldsc`, and `mixer`. No aliases or alternative spellings are accepted. When
omitted, formats come from YAML.
`--variant-id-type rsid` extracts an rsID from the VCF `ID` field;
`--variant-id-type unique` constructs the configured chromosome-position-REF-ALT
identifier. `variant_identifiers.target_types` can set different conventions
for different outputs. MAGMA pipeline mode inspects BIM field 2 and sets only
the MAGMA target automatically.

Duplicate-ID handling is shared by every built-in format and the custom table.
`--duplicate-id-policy` accepts `exclude_all`, `error`, `most_significant`,
`highest_maf`, or `highest_info`. Argparse owns no default: when the option is
omitted, each target uses `variant_identifiers.target_duplicate_policies` and
then the canonical `variant_identifiers.default_duplicate_policy`. The packaged
default is `exclude_all` for every target.

Exact repeated records are collapsed first. Supported formatter-side reference
matching—currently LDSC `--merge-alleles`—runs next when supplied. Other modules
retain their own downstream reference reconciliation. The resolved policy then
handles only conflicting duplicate-ID groups that remain. `exclude_all` removes the
complete group, and `error` stops. Ranked policies retain a row only when it has
one strictly greatest valid ranking value: largest configured `-log10(P)`,
largest `min(EAF, 1-EAF)`, or largest configured INFO value. A tied maximum or
missing ranking value excludes the complete group; the schema-validated
`duplicate_rank_tolerance` treats numerically equivalent floating-point ranks
as ties, and input order is never a fallback. `most_significant` is an explicit
user choice rather than the default
because selecting association results by P value can introduce ascertainment
bias. Every target records the policy and its exact resolution counts.

`--merge-alleles PATH` is an optional LDSC-formatting reference, normally
`w_hm3.snplist`. When supplied, the formatter retains only records whose rsID
and strand-unambiguous allele pair match one reference row. Allele order and
strand complement are recognized. This resolves a duplicated rsID when exactly
one input record is allele-compatible. Zero matches are excluded; multiple
compatible records are passed to the shared configured duplicate policy.
Without this option, direct LDSC formatting applies the same shared policy
directly. It never keeps the first duplicated row arbitrarily. The
heritability pipeline requires this option and passes the same file to both the
formatter and `munge_sumstats.py`.

`--custom-output FILE` activates the additional custom table. `--id NAME` is
mandatory and uses the same `--variant-id-type rsid|unique` policy. Each other
option requests one field and supplies its output header:

| Option | Value written |
|---|---|
| `--id NAME` | Selected rsID or configured chromosome-position-REF-ALT unique ID |
| `--chr NAME` | Canonical chromosome from the harmonised VCF |
| `--pos NAME` | One-based position |
| `--ref NAME` | Other/non-effect allele |
| `--alt NAME` | Effect allele |
| `--beta NAME` | Effect estimate |
| `--se NAME` | Standard error |
| `--z NAME` | Z statistic |
| `--lp NAME` | Input `-log10(P)` |
| `--p NAME` | Raw P converted from `LP`, bounded by the configured minimum |
| `--eaf NAME` | Effect-allele frequency |
| `--maf NAME` | `min(EAF, 1-EAF)` |
| `--n NAME` | Total sample size |
| `--neff NAME` | Effective sample size |
| `--n-case NAME` | Case sample size |
| `--n-control NAME` | Control sample size |
| `--info NAME` | Imputation INFO score |

There is no permissive missing-value option. Every requested field is required
for every written row; rows with missing, non-finite, or scientifically invalid
requested values are excluded and counted. Output header names must be unique.

Resume reuses outputs only after provenance, content, configuration, and
freshness validation. If a completed run directory is copied, resume validates
the copied files at the currently configured destinations, rebases every
returned artifact and manifest path to that directory, and rejects a recorded
artifact whose relative filename does not match the current configuration.
The detailed formatter HTML report is a PostGWAS-owned, SHA-256-fingerprinted
artifact under the same policy. Its links to formatter outputs are relative to
the configured output root, so they remain valid after a completed pipeline
stage is published or a validated run directory is copied.
When LDSC reference selection is active, the completion manifest also validates
the merge-alleles path, size, and SHA-256 digest.
If every formatter artifact declared by a matching manifest is absent, formatter
logs the stale completion record and regenerates all selected formats. If only
some artifacts are missing, or any surviving artifact has changed, formatter
stops and requires explicit `--overwrite` after review.
`--overwrite` takes precedence. Output schemas,
filenames, chromosomes, minimum representable P value, and MiXeR QC come from
the generated [Configuration Defaults](../reference/configuration-defaults.md).

## Processing steps

The formatter validates selection, renders every selected output destination,
and rejects filename collisions before VCF extraction. This preflight includes
named outputs and every configured chromosome partition. It then optionally
validates a completion manifest, validates the PostGWAS header contract,
extracts `N_ALT=1` records once, parses numeric
fields, selects each target's configured ID convention, and runs each exporter
in canonical order. Targets sharing the same resolved identifier type and
duplicate policy reuse one validated selection; the bounded cache retains only
selections with another consumer and releases them after the last consumer.
Every final candidate is checked for missing, empty, or duplicated selected IDs
before its exporter can write a scientific artifact. Before an LDSC export with
`--merge-alleles`, the formatter instead joins selected rsIDs to the configured
reference columns, applies the strand-unambiguous allele check, and validates
that separate result without using the shared cache. The custom exporter reuses
the same one-pass extraction and does not trigger study-design inference. Study
design is inferred only when LDSC or MiXeR is selected,
because only those formatter contracts interpret sample size differently for
binary and quantitative traits. MAGMA, SuSiE, and FINEMAP validate only their
own configured fields and do not require case/control columns. Rows missing or
violating a target's required fields are excluded for that target and counted.
Before PRED-LD field validation, rows whose canonical chromosome is absent
from the configured `chromosomes` list are excluded and counted by chromosome
label.

During a fresh run, the terminal shows two kinds of measurable formatter work:
one stage validates and reads the harmonised GWAS-VCF, followed by one stage for
each requested downstream format. The first completed-stage outcome reports the
total biallelic input-variant count, validated genome build, and dataset/sample
identity embedded in the VCF. Each output stage reports the same input total,
the number written, and the number excluded. The completion summary then names
the exact prepared files, the required-value exclusions, any configured
minimum-P bounding, and the validated VCF provenance. These are formatter
retention counts only: a downstream tool can retain fewer variants after its
own LD-reference, analysis-scope, or mapping checks.

When a BIM reference is provided for downstream identifier selection, formatter
reports the identifier convention detected in BIM field 2 and the number of BIM
records inspected. This inspection determines whether the prepared identifier
column uses rsIDs or the configured chromosome-position-allele convention. It
does **not** claim that the GWAS variants overlap the BIM reference; the
downstream module applies and reports its configured reference-matching policy.

For a MAGMA pipeline, the subsequent MAGMA runner receives the exact filtered
`snp_loc_file` and `pval_file` returned by this formatter step. MAGMA preparation
validates them again and the external command uses `duplicate=error`; therefore
formatter-produced duplicates cannot be silently reconsidered later. The
MAGMA module's lowest-P consolidation applies only to independently prepared
direct-module inputs that repeat one coordinate/allele-consistent variant.

For a binary LDSC export, each valid written variant contributes
`N_CASE / (N_CASE + N_CONTROL)`. The formatter reduces those per-variant case
fractions using `ldsc_sample_prevalence.aggregation`: `median` is the default
and `mean` is the only alternative. `sum` is rejected because variant rows are
not independent participant groups. The result, method, number of variants,
and observed minimum/maximum are recorded in the completion manifest and
canonical log. The same value, aggregation formula, contributing variant count,
observed prevalence range, `N_CASE` range, and `N_CONTROL` range are printed in
the terminal completion summary for both a fresh export and a validated resume.
All three ranges use the exact valid variants written to the LDSC table rather
than all input VCF records. Quantitative traits explicitly report that sample
prevalence is not applicable. A completion manifest created before count ranges
were recorded continues to resume safely and reports that a one-time
`--overwrite` rerun is required to populate the missing range metadata.

### Runtime column report

Every formatter execution, including a validated resume, prints one row per
selected output. The report is generated from the resolved YAML schema and
shows:

- the exact resolved variant-ID type for that target (`rsID` or `unique`, where
  `unique` is the configured coordinate-and-allele identifier);
- every canonical source column and its saved output header (`source → saved`);
- whether P values are absent, retained as `-log10(P)`, or converted to raw
  `P = 10^-LP`;
- the exact frequency source and saved header, with either EAF retained
  unchanged or `MAF = min(EAF, 1-EAF)` applied; and
- every exact GWAS-VCF sample-size source and saved header, its total/effective/
  case/control meaning, any trait-specific formula, and how the downstream tool
  uses it. `copied` means no numerical transformation occurs in the formatter.

The same structured report is written to the canonical formatter log under
`formatter_saved_schema`. Therefore custom YAML column names, CLI custom-column
names, and binary or quantitative LDSC sample-size columns are reported as
actually resolved for that run rather than described using fixed defaults.

### Terminal completion summary

The opening plan states that formatter creates downstream input files and does
not run the downstream analysis. The completion summary uses four visible
levels: the formatter result, major sections, target-specific subsections, and
their aligned fields. Each selected target contains separate variant-accounting
and identifier-handling subsections; the LDSC target also contains the GWAS-VCF
case-fraction subsection.
Reports and logs form a separate major section. Counts that represent a subset
include a percentage and state their denominator. For LDSC reference selection,
the summary reports the total number of unique rsIDs in the supplied LDSC
SNP/allele reference, the fraction of
GWAS-VCF variants retained, the fraction of reference variants represented in
the formatter output, absent reference rsIDs, allele mismatches, and duplicate
groups resolved by reference matching.

The label `Variants written to LDSC files` describes variants written to the
formatter's LDSC input table after its configured rsID, allele, duplicate, and
required-field checks. It does not claim that every variant will enter the
heritability regression: the later LDSC `munge_sumstats.py` stage may apply
additional INFO, MAF, sample-size, and statistical-quality filters. The supplied
LDSC SNP/allele reference is also kept distinct from the chromosome-split
LD-score reference and weight files used later by `ldsc.py`.

### Detailed HTML report

Every successful fresh formatter run writes the self-contained report configured
by `runtime.html_report_file`. The packaged path is
`reports/{dataset_id}_formatter_report.html`. A validated resume reuses the exact
checksum-matched report rather than regenerating it from potentially different
presentation code.

The report is rendered only from the same validated result metadata, resolved
schema, and paths used by the terminal, canonical log, and completion manifest;
it does not reopen the GWAS-VCF, exported tables, or LDSC reference. It contains:

- GWAS-VCF and sample-count completeness evidence used for trait inference;
- a target overview with retained/excluded counts and percentages;
- disjoint exclusion accounting for every selected formatter target;
- complete selected-identifier, reference-matching, duplicate-resolution, and
  final uniqueness evidence;
- the total supplied LDSC SNP/allele reference size, matched coverage, absent
  rsIDs, allele mismatches, and reference-resolved duplicate groups;
- the GWAS-VCF case-fraction formula, aggregation, contributing variants, value
  range, and case/control count ranges;
- every configured canonical-source to saved-column mapping and transformation;
- P-value, allele-frequency, and sample-size semantics for each downstream tool;
- the resolved formatter policies and canonical GWAS-VCF input contract; and
- relative links to every generated artifact plus the canonical log, resolved
  configuration, checksum completion manifest, input VCF, and executable
  provenance.

With the default configuration, the identifier and statistical representations
are:

| Output | Variant ID type | P value saved | Frequency saved | Sample size saved |
|---|---|---|---|---|
| MAGMA | `rsID` | `LP → P` as raw P | Not saved | Total N: `FORMAT/SS → N_COL`, copied per variant |
| GCTA gene | `rsID` | `LP → P` as raw P | `EAF → freq` unchanged | Total N: `FORMAT/SS → N`, copied per variant |
| SuSiE | `rsID` | `LP → LP` unchanged | Not saved | `FORMAT/NEF → NEF`; the locus median is used as SuSiE `n` |
| FINEMAP | `rsID` | Not saved | `EAF → maf` using `min(EAF, 1-EAF)` | `FORMAT/NEF → NEF`; the rounded locus median is used as `n_samples` |
| PRED-LD | `rsID` | `LP → LP` unchanged | `EAF → AF` unchanged | `FORMAT/NCO → NC` and `FORMAT/SS → SS`; carried for re-harmonisation, not consumed by PRED-LD |
| LDSC binary | `rsID` | `LP → P` as raw P | `EAF → FRQ` unchanged | Cases: `FORMAT/NC → N_CAS`; controls: `FORMAT/NCO → N_CON` |
| LDSC quantitative | `rsID` | `LP → P` as raw P | `EAF → FRQ` unchanged | Total N: `FORMAT/NCO → N` |
| MiXeR binary | `rsID` | Not saved | Not saved | Effective N: `FORMAT/NEF → N`, where `NEF = 4/(1/NC + 1/NCO)` |
| MiXeR quantitative | `rsID` | Not saved | Not saved | Total N: `FORMAT/NEF → N`, where `NEF = NCO` |
| Custom CLI table | `rsID` | Exact `--p`/`--lp` source and requested header | Exact `--eaf`/`--maf` source and requested header | Exact requested sample-size source and header |

## Outputs

- MAGMA: `<dataset>_magma_snp_loc.tsv` and `<dataset>_magma_p_values.tsv`.
- GCTA gene: `<dataset>_gcta.ma`, shared by fastBAT and mBAT-combo.
- SuSiE: `<dataset>_susie.tsv`.
- FINEMAP: `<dataset>_finemap.tsv`.
- PRED-LD: `pred_ld/<dataset>_chr<chromosome>_pred_ld_input.tsv` for configured chromosomes.
- LDSC: `<dataset>_ldsc_input.tsv`.
- MiXeR: `<dataset>_mixer.sumstats.gz`.
- Custom: the relative filename supplied to `--custom-output`.
- Detailed HTML: `reports/<dataset>_formatter_report.html` by default; the path
  comes from `runtime.html_report_file` in canonical YAML.
- Provenance: `logs/<dataset>_formatter.log`,
  `run_metadata/resolved_config.yaml`, and
  `run_metadata/formatter_completion.yaml`.

## QC and logs

For each target, review the selected ID type, identifier exclusions, rows
in/out, exclusions for missing or invalid values,
P values bounded at the configured numeric minimum, written schema, inferred
study design when required, sample-size mode, and output fingerprints. PRED-LD additionally
records `rows_excluded_unconfigured_chromosome` in its result and lists each
excluded canonical chromosome with its row count in the canonical log.
Every format records duplicate selected-ID groups, the resolved policy, and the
number of exact, ranked, ambiguous, and excluded rows. The terminal completion
summary reports detected duplicate groups and rows, reference-resolved groups,
groups not retained by reference matching, groups still ambiguous after
reference matching, exact rows collapsed, policy-resolved groups, rows removed
by policy, final retained rows, percentages, and the successful uniqueness
invariant. The canonical log and completion metadata also
record whether a target reused an already validated identifier selection. LDSC
reference selection additionally records unmatched rsIDs, allele mismatches,
groups resolved uniquely by alleles, and groups still ambiguous after reference
matching.
The custom result additionally records its ordered field roles, output headers,
identifier convention, and missing/invalid requested-field exclusions.

## Interpretation

Allele direction is preserved: GWAS-VCF ALT is the effect allele. FINEMAP gets
minor allele frequency; MAGMA and GCTA mBAT receive total sample size. SuSiE and
FINEMAP use locus-median `NEF`, which is effective N for binary traits and total
N for quantitative traits. LDSC uses case/control counts for binary traits or N
from NCO for quantitative traits; MiXeR likewise uses binary effective N or
quantitative total N from `NEF`.

## Common problems

No format selected, colliding resolved output filenames, no usable IDs of the
selected type, duplicated IDs for a target requiring strict uniqueness, every
LDSC rsID belonging to a duplicate group, missing VCF fields, incomplete
per-variant case/control counts, stale outputs, a changed resolved config,
manifest artifact paths that do not match the current configured filenames, or
zero rows satisfying one target's scientific contract. With
`--merge-alleles`, duplicated reference rsIDs and invalid or strand-ambiguous
reference allele pairs are also fatal. Multiple compatible VCF records for one
reference rsID are handled by the configured duplicate-ID policy.

## Limitations

Formatting does not run the external analyses and cannot make an inappropriate
reference panel compatible. A custom table guarantees the requested formatter
semantics but does not assert that the table satisfies an arbitrary external
tool's complete input contract. When LDSC or MiXeR requires study design, it is
inferred from count values rather than a separate VCF trait declaration.

## Scientific references

See [Scientific References](../reference/scientific-references.md) and the
documentation for the downstream tool receiving each export.
