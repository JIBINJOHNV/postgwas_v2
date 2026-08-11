# Harmonisation Sample Sheet

The sample sheet is a CSV or TSV with one dataset per row. Its header is
normalized for case, whitespace, Unicode compatibility, and hyphens; unknown
fields are rejected. Referenced data files are checked for existence and size
during sample-sheet preflight but opened only by later execution preflight.

## Required identity and file fields

| Field | Requirement |
|---|---|
| `config_version` | Must be `2`. |
| `dataset_id` | Unique; letters, numbers, `.`, `_`, and `-`, beginning with a letter or number. |
| `input_file` | Existing, non-empty regular file. Relative paths are resolved from the sample-sheet directory. |
| `trait_type` | `quantitative`, `case_control`, or `auto`. |
| `delimiter` | `tab`, `comma`, `semicolon`, `space`, `whitespace`, or `auto`. |

## Coordinate and allele mappings

Provide either `chromosome_column` plus `position_column`, or
`chromosome_position_column`. `effect_allele_column` and
`other_allele_column` are required. `variant_id_column` is optional because an
identifier can be constructed after validation.

Allele order is scientifically meaningful: the effect estimate and effect
allele frequency refer to `effect_allele_column`.

## Effect and P-value mappings

Provide `effect_column` or `z_score_column`. `effect_type` accepts `beta`,
`odds_ratio`, or `auto`; odds ratios are converted to the log scale during
harmonisation. `standard_error_column` is optional because PostGWAS can derive a
standard error from valid statistics under configured policies.

`p_value_column` is required. Declare `p_value_type` as `raw`, `neglog10`,
`negln`, or `auto` so transformed P values are not treated as raw values.

## EAF and INFO mappings

Choose exactly one source for study EAF. INFO uses the explicit priority shown
below:

| Concept | Internal source | External source | Resolution rule |
|---|---|---|---|
| Study effect-allele frequency | `effect_allele_frequency_column` | `external_eaf_file` plus `external_eaf_column` | Exactly one source is required. |
| Study imputation quality | `imputation_info_column` | `external_info_file` plus `external_info_column` | Internal INFO takes priority when both are listed; otherwise use external INFO. With neither, explicitly supply `--fixed-info VALUE` or preflight fails. |

An external path without its column name is incomplete. The comparison AF
panel configured for harmonisation is not a substitute for either study source.
An ignored external INFO source is reported rather than silently combined with
the internal values.

An external EAF or INFO source may be either one existing table containing all
chromosomes or an explicit per-chromosome path template. Use `{chromosome}` for
the chromosome label and optionally `{build}` for the inferred build, for
example `/references/GRCh37_panel_chr{chromosome}.tsv.gz`. Do not enter only the
common filename prefix: PostGWAS does not guess which matching file belongs to
each chromosome.

## Sample-size mappings

For a quantitative trait, provide `control_count_column` or `control_count`;
this field represents total analyzed N. Do not provide a case count.

For a case-control trait, provide controls using `control_count_column` or
`control_count`, and cases using `case_count_column` or `case_count`. Counts
must be positive whole numbers. Effective N is
`4 / (1/Ncase + 1/Ncontrol)`.

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
python -m postgwas.modules.harmonisation.sample_sheet_generator \
  --input-directory /path/to/summary-statistics \
  --output studies.csv
```

The generator reads only enough records to determine the delimiter and header.
Every supported input file becomes one row in the same output CSV, ordered by
filename. Selected column mappings preserve the source header spelling exactly,
including a leading `#` in names such as `#chr` or `#chrom`. The terminal summary
reports the number of ready rows and lists only datasets requiring attention.
The complete recognized file suffix is removed from `dataset_id`; for example,
`study.vcf.gz` becomes `study`, not `study.vcf`.
An unresolved EAF mapping is written as `NA` so the user can add either the
exact internal frequency column or an external EAF file/column pair. Unresolved
chromosome, position, allele, effect, or p-value mappings still stop generation
without publishing a partial multi-study sheet.

It writes the same fields as the validated sample-sheet model and sets
`trait_type`, `effect_type`, and `p_value_type` to `auto` for every dataset.
Header aliases select sample-size, effect, and p-value columns only; they do not
establish the study design or statistical scales. Dataset-level harmonisation
makes those decisions from the complete data. Edit these fields manually when
the study documentation provides an authoritative declaration.

Missing EAF and missing sample size are written as `NA` in the draft. The
terminal summary identifies the affected dataset and explains how to complete
its EAF source. The normal harmonisation sample-sheet validator remains strict:
it rejects the draft until exactly one internal or external EAF source is
provided, and it never invents a fixed frequency.

Missing sample size is written as `NA` under the default
`sample_sheet_generator.on_missing_sample_size: write_draft` policy. The
terminal warning identifies the exact fields that must be completed, and the
normal harmonisation sample-sheet validator still rejects the draft until they
are supplied. Set the policy to `fail` to retain strict generation. Other
unresolved required mappings stop generation without publishing partial
output, except for EAF as described above.

MAF-like and allele-ambiguous frequency headers are retained with visible
review warnings; reference-panel frequency columns are never relabelled as
study EAF. An ALT-labelled frequency is accepted as EAF only when the selected
effect-allele column is also ALT-labelled. Review every generated allele,
frequency, effect, and sample-size mapping before starting harmonisation. Add
external EAF/INFO sources manually when the measurements are not present in the
study table.

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
missing alternative mappings, invalid trait/count combinations, extra columns,
or missing/empty referenced files. Validation errors identify the CSV row,
dataset and exact fields that must be corrected without exposing internal
validation-library output. A header can pass this stage and still fail
execution if a mapped column is absent from the actual summary-statistics file.
