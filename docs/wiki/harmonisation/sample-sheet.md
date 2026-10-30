# Harmonisation Sample Sheet

The sample sheet is a CSV or TSV with one dataset per row. Its header is
normalized for case, whitespace, Unicode compatibility, and hyphens. Extra
columns, such as `phenotype`, `pubmed_id`, or custom review notes, are allowed
and ignored during harmonisation. They remain in your source sheet and are
not included in the normalized analysis row. Required fields and the values
of recognized fields are still validated. Referenced data files are checked for existence and size
during sample-sheet preflight but opened only by later execution preflight.

## Required identity and file fields

| Field | Requirement |
|---|---|
| `config_version` | Must be `2`. |
| `dataset_id` | Unique; letters, numbers, `.`, `_`, `+`, and `-`, beginning with a letter or number. |
| `input_file` | Existing, non-empty regular file. Relative paths are resolved from the sample-sheet directory. |
| `trait_type` | `quantitative`, `case_control`, or `auto`. |
| `delimiter` | `tab`, `comma`, `semicolon`, `space`, `whitespace`, or `auto`. |

## Coordinate and allele mappings

Provide either `chromosome_column` plus `position_column`, or
`chromosome_position_column`. `effect_allele_column` and
`other_allele_column` are required. `variant_id_column` is optional because an
identifier can be constructed after validation.

Positions are ordinary 1-based GWAS coordinates and must resolve to finite
whole numbers of at least 1. Integer-valued representations such as `100`,
`100.0`, and `1e2` are accepted without changing their value. Fractional,
non-finite, zero, negative, and unreadable positions are rejected as
`invalid_position`; PostGWAS never truncates a fractional coordinate.

Allele order is scientifically meaningful: the effect estimate and effect
allele frequency refer to `effect_allele_column`.

Before running a generated or copied sheet, compare its declarations with the
study's data dictionary: trait design, effect type and scale, P-value scale,
allele order, frequency meaning and sample-size definition. Set known
`trait_type`, `effect_type` and `p_value_type` explicitly; retain `auto` only
when inference is actually needed. An automatically recognized column name is
not evidence that its values have the intended meaning.

## Effect and P-value mappings

Provide `effect_column` or `z_score_column`. `effect_type` accepts `beta`,
`odds_ratio`, or `auto`; odds ratios are converted to the log scale during
harmonisation. `standard_error_column` is optional because PostGWAS can derive a
standard error from valid statistics under configured policies.

`p_value_column` is required. Declare `p_value_type` as `raw`, `neglog10`, or
`auto` so transformed P values are not treated as raw values. Convert signed
logarithms and natural-log representations to raw P or non-negative -log10 P
before harmonisation.

## EAF and INFO mappings

Choose exactly one source for study EAF. INFO uses the explicit priority shown
below:

| Concept | Internal source | External source | Resolution rule |
|---|---|---|---|
| Study effect-allele frequency | `effect_allele_frequency_column` | `external_eaf_file` plus `external_eaf_column` | Exactly one source is required. |
| INFO value and provenance | `imputation_info_column` | `external_info_file` plus `external_info_column` | Internal INFO takes priority when both are listed; otherwise use external INFO. With neither, explicitly supply `--fixed-info VALUE` or preflight fails. |

A selected external source needs both its path and column name. The comparison
AF panel configured for harmonisation is not a substitute for either study
source. An ignored external INFO source is reported rather than silently
combined with the internal values. External INFO is a reference proxy, and
fixed INFO is a user-assigned value; neither establishes the study's measured
imputation quality. Read their provenance before applying an INFO filter.

An external EAF or INFO source may be either one existing table containing all
chromosomes or an explicit per-chromosome path template. Use `{chromosome}` for
the chromosome label and optionally `{build}` for the inferred build, for
example `/references/GRCh37_panel_chr{chromosome}.tsv.gz`. Do not enter only the
common filename prefix: PostGWAS does not guess which matching file belongs to
each chromosome.

`external_eaf_column` declares an ALT/effect-allele frequency, not MAF. Before
allele-swapped matches apply `1-AF`, PostGWAS screens the deduplicated raw column.
A non-MAF-like distribution retains that declaration. A MAF-like distribution
must be confirmed against the independently configured
`modules.harmonisation.default_eaf` panel. Confirmed MAF or inconclusive evidence
stops; PostGWAS does not auto-convert it. If both sources resolve to the same
physical file, choose a different `modules.harmonisation.default_eaf.source` so
the comparison is independent. Before deciding EAF versus MAF, PostGWAS also
requires the supplied and default panels to have compatible folded-MAF spectra
under the YAML-configured correlation and mean-absolute-difference limits. This
checks whether the chosen default population panel is suitable; it does not use
correlation alone to infer allele direction.

## Sample-size mappings

For a quantitative trait, provide `control_count_column` or `control_count`;
this field represents total analyzed N. Do not provide a case count.

For a case-control trait, provide controls using `control_count_column` or
`control_count`, and cases using `case_count_column` or `case_count`. Counts
must be positive whole numbers. Effective N is
`4 / (1/Ncase + 1/Ncontrol)`.

### A case-control file containing only precomputed Neff

The version-2 sample sheet has no separate effective-N field. If the source
contains a valid, already calculated Neff but does not contain case and control
counts, mapping that Neff column to `control_count_column` preserves the number
for downstream consumers that require effective N only. MiXeR is such a
consumer: its
[official input guidance](https://github.com/precimed/mixer)
defines case-control N as effective N.

Use this sample-sheet mapping only when the source column is explicitly defined
as precomputed effective N by the study documentation:

```csv
trait_type,control_count_column,case_count_column
auto,NEFF,NA
```

Do not declare the binary source quantitative. An explicit `case_control`
declaration requires case counts, so `auto` is the only current compatibility
route when those counts are unavailable.

This mapping is a compatibility workaround, not a complete representation of
the study design. Without a case column, downstream formatting classifies the
study as quantitative. The resulting VCF records the supplied value as total
sample size and control count as well as effective sample size. Consequently,
do not use the workaround when an analysis needs any of the following:

- actual total sample size;
- separate case and control counts;
- sample prevalence or automatic case-control classification;
- liability-scale LDSC without separately supplied sample prevalence.

The input column must already be Neff. A binary study's total N alone, or its
actual control count alone, is not Neff and cannot determine Neff without the
case/control ratio. When those downstream meanings matter, obtain and map both
case and control counts instead.

## Inference and aliases

Missing `trait_type`, `effect_type`, `p_value_type`, and `delimiter` values
normalize to `auto`. Common aliases such as `binary`, `OR`, `pvalue`, `-log10p`,
and `csv` are normalized. An invalid inferable value becomes `auto` with a
warning; later scientific preflight decides whether the data provide enough
evidence to proceed.

## Generate a draft from GWAS headers

The configured header aliases can generate a version-2 sample-sheet draft for
every supported summary-statistics file in one directory:

```console
postgwas sample_sheet_generator \
  --input-directory /path/to/summary-statistics \
  --output studies.csv
```

### Optional external sources for missing mappings

These options are opt-in; the basic generator command needs neither pair.

| What you want to supply | Optional CLI pair |
|---|---|
| External EAF only | `--external-eaf-file PATH --external-eaf-column NAME` |
| External INFO only | `--external-info-file PATH --external-info-column NAME` |
| Both | Supply both pairs shown below. |

Each chosen pair requires both values after CLI/YAML resolution. Replace the
column names with the exact headers in your reference tables. Without a
configured external pair, missing source fields remain `NA` and are reported
for attention; the generator does not download or invent a replacement.

```bash
postgwas sample_sheet_generator \
  --input-directory /path/to/summary-statistics \
  --output studies.csv \
  --external-eaf-file '/references/{build}/eaf_chr{chromosome}.tsv.gz' \
  --external-eaf-column AF \
  --external-info-file '/references/{build}/info_chr{chromosome}.tsv.gz' \
  --external-info-column INFO
```

Each optional file/column pair applies to every retained study lacking that
internal mapping. Either pair can be supplied alone. Internal frequency
(including MAF/reference-frequency candidates) and INFO mappings are never
replaced. A mapped column with missing or invalid values is still an internal
mapping; value validation occurs during harmonisation.

The four corresponding settings are
`modules.harmonisation.sample_sheet_generator.external_eaf_file`,
`external_eaf_column`, `external_info_file`, and `external_info_column`, all
unset by default. CLI values override YAML before pair validation. Relative
CLI paths resolve from the working directory; relative YAML paths resolve
from the supplied run-config directory. Absolute paths, retaining template
placeholders, are written to the CSV.

Generation validates pairs, template syntax, and matching non-empty readable
files once, and excludes those reference files from study discovery. It does
not read the complete external tables or infer study builds. Harmonisation
preflight checks the selected build/chromosome files and configured external
column mappings. `{build}` and `{chromosome}` select files; they do not lift
over coordinates or search alternative builds after a failed match.

The review TSV adds `eaf_source` and `info_source` (`internal`, `external`, or
`missing`). Missing-source flags account for supplied external pairs. The
terminal separately counts external assignments and rows without internal
INFO. Assigned external sources receive review notes: their build, population,
and provenance need confirmation for each study, and reference INFO is not
study-measured imputation quality. The v2 sample-sheet columns are unchanged.

The generator reads only enough records to determine the delimiter and header.
Supported inputs include `.tbl`, `.tbl.gz`, `.tbl.zip`, and other supported
text tables ending in `.gz` or `.zip`. A ZIP archive must contain exactly one
data file (directories and `__MACOSX/` entries are ignored). The same header
screening applies regardless of extension.
Every valid supported input file becomes one row in the same output CSV,
ordered by filename. Selected column mappings preserve the source header
spelling exactly, including a leading `#` in names such as `#chr` or `#chrom`.
The terminal prints aggregate counts only: scanned/written/rejected files,
ready/attention rows, missing INFO/sample size/frequency/mappings, mapped
frequency/coordinates, and rows with both case/control columns. It lists report
paths instead of individual datasets or warnings. Counts describe retained
header mappings; they do not establish study design or frequency type.
The complete recognized file suffix is removed from
`dataset_id`; for example, `study.vcf.gz` becomes `study`, not `study.vcf`.
The filename-derived ID must satisfy the same rules as a manually entered ID.
IDs containing `+`, such as `EVER_SMOKER_GWAS_MA_UKB+TAG`, are preserved.
An invalid ID, such as a filename stem containing spaces, follows the
configured per-file rejection policy instead of aborting CSV publication with
a validation traceback. The generator does not rename source files or silently
rewrite IDs. To include a rejected study manually, use a valid `dataset_id`
while retaining its original path in `input_file`.
The generator follows the draft-first approach of the
[original mapping script](https://github.com/JIBINJOHNV/postgwas/blob/main/src/postgwas/scripts/create_sumstat_map_pl.py)
while retaining the version-2 columns. It chooses the first configured alias
for each field, preserves the exact input header, and records alternatives for
review. Missing mappings are written as `NA` instead of rejecting an otherwise
recognizable study. It does not create harmonisation result directories.

The canonical `sample_sheet_generator.screening_columns` lists six core roles:
chromosome, position, effect allele, other allele, effect, and p-value.
`minimum_mapped_columns: 3` retains files with at least three mapped roles,
matching the original script's inclusion threshold. Combined coordinates satisfy
both coordinate roles, and a Z-score can satisfy the effect role. This is a
file-screening heuristic, not scientific validation. The normal harmonisation
loader still requires complete mappings before analysis.

If normal delimiter detection fails, the generator can rescue mixed tabs and
spaces by writing a separate, validated copy. For example,
`ILAE3_CAE_final.tbl` becomes `prepared_inputs/ILAE3_CAE_final.tbl.converted.txt`
beside the output CSV. The CSV records its **absolute path**, keeps the original
dataset ID, and uses `delimiter: space`. Originals are never modified.

The fallback applies to all discovered supported extensions, including
compressed inputs, which are streamed through decompression first. It performs
the equivalent of `tr '\t' ' ' | tr -s ' '` and validates every data row.
Quoted text, positional empty tab fields, boundary whitespace, invalid text,
and inconsistent converted field counts are rejected. Files failing only
column-name mapping or other non-delimiter checks are not converted. Normal
files retain the existing sampled delimiter detection; this is not a full-file
structural audit of every input.

The terminal reports conversions attempted, accepted and rejected on one line.
The review TSV retains the original `input_file` and adds `prepared_input_file`,
`conversion_status` and `converted_rows`, with transformation details in
`review_notes`. Converted rows remain flagged for review. Identical existing
copies are reused; differing files or symlink destinations are refused.
Generated copies are excluded from subsequent study discovery.

Configure `modules.harmonisation.sample_sheet_generator` with
`convert_mixed_whitespace: true` (default), `prepared_input_directory: prepared_inputs`
and `converted_file_suffix: .converted.txt`. The directory is a single folder
name relative to the CSV output directory. Set `convert_mixed_whitespace: false`
to disable rescue. Conversion fixes separators only; it does not add missing
frequency, sample-size or INFO mappings.

Two separate TSV reports accompany the CSV:

- `<output>.mapping_review.tsv`: one row per candidate file, with `dataset_id`,
  `input_file`, `status` (`ready`, `needs_attention`, or `rejected`),
  `missing_mappings`, `missing_info`, `missing_sample_size`, `missing_frequency`,
  and `review_notes`. Missing-field flags are `True` or `False` for retained
  rows; rejected rows leave unassessed flags blank and include the rejection
  reason. Notes include MAF uncertainty, reference-frequency provenance,
  ALT alignment, and alternative matching columns.
- `<output>.rejected_files.tsv`: full paths and reasons for files that cannot
  be read or fail the configured screening/structural checks.

The report suffixes are configured by `review_report_suffix` and
`rejection_report_suffix`. Reports and the output sheet are excluded from
subsequent input discovery. `ready` means no header-mapping issue was detected;
it does not establish scientific validity. Generation fails if no candidate
passes screening, but still writes both diagnostic reports.
`failure_policy: fail_all` requires every candidate to pass before a sample
sheet is written; it still writes the diagnostic reports on failure;
the default `write_valid` writes the retained rows and reports rejected files.

`trait_type`, `effect_type`, and `p_value_type` remain `auto`. Header matching
does not establish study design or statistical scale. Edit these fields when
study documentation supplies an authoritative declaration.

Frequency candidates including EAF, MAF, ambiguous frequency, ALT frequency,
and `EAF_HRC` are retained in `effect_allele_frequency_column`. Explicit EAF
aliases take priority, followed by ALT, MAF, ambiguous, and reference-frequency
aliases; the configured order within each category breaks ties. Prefix matches
are the fallback. MAF and reference-frequency candidates receive review notes,
not generation errors. Retaining a header does not establish study provenance
or allele orientation; harmonisation's existing frequency/reference checks
remain unchanged and can still stop on inconclusive evidence.

Missing sample size is written as `NA` under `on_missing_sample_size: write_draft`.
Set this to `fail` to reject that file during generation. Explicit controls take
priority over total/effective N. If cases and total N are present without a
control-count column, that total/effective column is not used for controls.
If no control/sample-size mapping was selected, `Weight` (case-insensitive)
is mapped to `control_count_column` as a last-resort draft fallback, including
when a case-count column exists. Its original header spelling is retained.
The row remains flagged for attention: confirm whether the values represent
total N, effective N, or actual controls before harmonisation. This mapping
does not convert weights into counts or alter numerical sample-size validation.

The alias key is `sample_sheet_generator.column_aliases.fallback_sample_size_column`
(default `[WEIGHT]`), replacing `invalid_sample_size_column`; update old custom
YAML overrides to the new key. Explicit controls and eligible total/effective
N mappings take precedence. If none of these or the fallback is present,
controls remain unresolved. These choices and missing INFO are listed in the
review table. Complete the mappings or supply appropriate
external EAF/INFO sources or `--fixed-info` before harmonisation.

## Examples

Use the maintained templates:

```text
examples/configs/harmonisation/sample_sheet_quantitative.csv
examples/configs/harmonisation/sample_sheet_case_control.csv
```

These contain placeholder paths and column names. Copy one, then replace every
study-specific value.

## Preflight failures

The sheet is rejected for wrong version, duplicate dataset IDs, invalid IDs,
missing required fields or alternative mappings, invalid trait/count
combinations, or missing/empty referenced files. Extra metadata columns do not
cause rejection. Validation errors identify the CSV row,
dataset and exact fields that must be corrected without exposing internal
validation-library output. A header can pass this stage and still fail
execution if a mapped column is absent from the actual summary-statistics file.
