# Input Data

PostGWAS accepts two broad forms of input: raw summary-statistics tables for
harmonisation, and harmonised GWAS-VCF or module-specific derivatives for
downstream analysis.

## Raw summary statistics

The harmonisation sample sheet uses `config_version: 2` and one row per dataset.
Each row identifies the source file and maps its study-specific columns.

### Example summary-statistics file

PostGWAS does not require the source columns to use these exact names. The
sample sheet maps the names used by your study to the concepts required by
harmonisation. For example, a tab-separated quantitative-trait file could look
like this:

```text
CHR	BP	SNP	A1	A2	EAF	BETA	SE	P	INFO
1	10177	rs_example_1	A	G	0.425	0.012	0.006	0.0455	0.98
1	10352	rs_example_2	T	C	0.118	-0.031	0.010	0.00194	0.96
2	20001	rs_example_3	C	T	0.276	0.004	0.008	0.617	0.91
```

These rows are illustrative rather than real reference variants. In this
example:

- `CHR` and `BP` identify the genomic coordinate in the declared genome build.
- `A1` is the effect allele and `A2` is the other allele.
- `EAF` is the frequency of `A1`, not necessarily ALT-allele frequency.
- `BETA` and `SE` are defined with respect to `A1`; `P` tests that association.
- `INFO` is the study's imputation-quality value.

The corresponding minimal sample sheet can use an analyzed sample size that is
constant across variants:

```csv
config_version,dataset_id,input_file,trait_type,chromosome_column,position_column,variant_id_column,effect_allele_column,other_allele_column,effect_allele_frequency_column,effect_column,effect_type,standard_error_column,p_value_column,p_value_type,control_count,imputation_info_column,delimiter
2,height_gwas,summary/height_gwas.tsv.gz,quantitative,CHR,BP,SNP,A1,A2,EAF,BETA,beta,SE,P,raw,500000,INFO,tab
```

`summary/height_gwas.tsv.gz` is resolved relative to the directory containing
the sample sheet. For a per-variant sample size, replace `control_count` with
`control_count_column` and give the name of the N column in the input file.

### Sample-sheet columns

Use the sample-sheet column names shown below. Do not add unsupported columns.

| Sample-sheet column | Requirement | Meaning and accepted values |
|---|---|---|
| `config_version` | Required | Sample-sheet schema version. It must be `2`. |
| `dataset_id` | Required | Unique output/study identifier. It must begin with a letter or number and may contain letters, numbers, `.`, `_`, and `-`. |
| `input_file` | Required | Raw summary-statistics file. A relative path is resolved from the sample-sheet directory; the file must exist and be non-empty. |
| `trait_type` | Optional/inferred | `quantitative`, `case_control`, or `auto`. Aliases include `continuous` and `binary`. Use an explicit type whenever it is known. |
| `chromosome_column` | Conditionally required | Name of the input chromosome column. Supply it together with `position_column`, unless using `chromosome_position_column`. |
| `position_column` | Conditionally required | Name of the input base-pair-position column. Supply it together with `chromosome_column`, unless using `chromosome_position_column`. |
| `chromosome_position_column` | Alternative coordinate mapping | Name of one input column containing a combined chromosome and position. Use this instead of, or in addition to, the separate coordinate pair; at least one complete coordinate mapping is required. |
| `variant_id_column` | Optional | Name of the input variant-ID column. Harmonisation can construct an identifier after coordinate and allele validation when this is absent. |
| `effect_allele_column` | Required | Input column containing the allele to which the effect, Z score, P value, and effect-allele frequency refer. |
| `other_allele_column` | Required | Input column containing the non-effect/other allele. |
| `effect_allele_frequency_column` | Exactly one EAF source | Input column containing study effect-allele frequency. Do not also supply external EAF fields. |
| `external_eaf_file` | Exactly one EAF source | External study-EAF file. It must be paired with `external_eaf_column` and cannot be combined with `effect_allele_frequency_column`. |
| `external_eaf_column` | Required with external EAF file | Name of the EAF value column in `external_eaf_file`. |
| `effect_column` | Conditionally required | Input effect-estimate column. Supply this or `z_score_column`; both may be mapped when present. |
| `effect_type` | Optional/inferred | `beta`, `odds_ratio`, or `auto`; `OR`, `odds ratio`, and `odds-ratio` normalize to `odds_ratio`. Odds ratios are converted to log-odds effects during harmonisation. |
| `standard_error_column` | Optional | Input standard-error column. PostGWAS can derive SE from other valid mapped statistics when the selected policy permits it. |
| `z_score_column` | Alternative effect mapping | Input Z-score column. Supply this when `effect_column` is absent. |
| `p_value_column` | Required | Input association P-value column. Its representation is declared by `p_value_type`. |
| `p_value_type` | Optional/inferred | `raw`, `neglog10`, `negln`, or `auto`. Recognized aliases include `pvalue`, `-log10p`, `mlogp`, and `-ln`. |
| `control_count_column` | At least one control/N source is required | Input column containing total analyzed N for quantitative traits or control N for case-control traits. Use this for per-variant counts. |
| `control_count` | At least one control/N source is required | Positive whole-number constant used as total analyzed N for quantitative traits or control N for case-control traits. Use this when the count is the same for all variants. |
| `case_count_column` | Required for case-control traits unless constant supplied | Input column containing case N. Use this for per-variant counts. Quantitative traits must not provide case-count fields. |
| `case_count` | Required for case-control traits unless column supplied | Positive whole-number constant case N. Use this when the count is the same for all variants. Quantitative traits must not provide case-count fields. |
| `imputation_info_column` | Exactly one INFO source | Input column containing study imputation quality. Do not also supply external INFO fields. |
| `external_info_file` | Exactly one INFO source | External study-INFO file. It must be paired with `external_info_column` and cannot be combined with `imputation_info_column`. |
| `external_info_column` | Required with external INFO file | Name of the INFO value column in `external_info_file`. |
| `delimiter` | Optional/inferred | `tab`, `comma`, `semicolon`, `space`, `whitespace`, or `auto`. Aliases include `\\t`, `csv`, and `comma-separated`. |

For case-control traits, effective sample size is calculated as
`4 / (1/Ncase + 1/Ncontrol)`. The sample sheet must therefore provide both a
case source and a control source, using either mapped columns or positive
whole-number constants.

See [Harmonisation Sample Sheet](../harmonisation/sample-sheet.md) for validation
rules, quantitative and case-control templates, and common errors.

## Delimiters and missing values

The declared delimiter may be `auto`, tab, comma, semicolon, space, or generic
whitespace. Automatic detection should still be verified against the validation
summary, particularly for files containing mixed whitespace or quoted fields.
Do not encode scientific missingness as a plausible numeric sentinel.

## Harmonised GWAS-VCF

The GWAS-VCF carries coordinates, alleles, study statistics, and metadata in a
form understood by downstream modules. Keep its index beside it where one is
created. Do not hand-edit VCF fields or headers between stages.

## Module-specific inputs

Formatting creates tool-specific tables for consumers such as LDSC, MAGMA,
GCTA gene analysis, fine-mapping, PrediXcan-style LD imputation workflows, and
MiXeR. A table formatted for one consumer is not automatically valid for
another. Follow the producing and consuming module pages together.

See [Input and Output Contracts](../core/input-output-contracts.md) for the
metadata that must remain consistent across these boundaries.
