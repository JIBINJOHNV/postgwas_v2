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

Choose exactly one source for each concept:

| Concept | Internal source | External source |
|---|---|---|
| Study effect-allele frequency | `effect_allele_frequency_column` | `external_eaf_file` plus `external_eaf_column` |
| Study imputation quality | `imputation_info_column` | `external_info_file` plus `external_info_column` |

An external path without its column name is incomplete. The comparison AF
panel configured for harmonisation is not a substitute for either study source.

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
or missing/empty referenced files. A header can pass this stage and still fail
execution if a mapped column is absent from the actual summary-statistics file.
