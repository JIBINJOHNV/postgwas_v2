# MAGMAcovar

## Purpose

MAGMAcovar runs MAGMA gene-property analysis using a completed MAGMA
`.genes.raw` file and a continuous gene-level covariate table.

## What the analysis does

MAGMA converts gene-association p-values to gene Z-scores and regresses them on
each supplied property while accounting for gene-gene correlation and MAGMA's
internal corrections. PostGWAS validates the inputs, runs MAGMA in an isolated
staging directory, validates every published `COVAR` row, and records the exact
command and resolved configuration.

The packaged standalone default is MAGMA's two-sided marginal gene-property
test. Choose a one-sided direction only for a prespecified directional
hypothesis. FUMA-style tissue specificity uses a one-sided `greater` test for
each tissue while conditioning on the cross-tissue `Average` property.

## When to use it

Use MAGMAcovar after MAGMA gene association when testing whether a continuous
gene property, such as expression in a tissue, is associated with gene-level
GWAS signal. A FLAMES pipeline uses this module to create tissue-relevance
weights. The original FLAMES README uses marginal two-sided tests; the separate
FUMA integration uses the conditional one-sided model described below.

## Input requirements

- A MAGMA `.genes.raw` file containing `# VERSION` metadata and readable gene
  rows.
- A covariate table with a header, gene IDs in column one, at least one numeric
  property, unique rows and columns, and `NA` for missing values.
- Matching gene identifiers across both inputs and at least the configured
  minimum number of overlapping genes.
- Variable properties within the allowed missingness threshold. PostGWAS rejects
  properties MAGMA would otherwise silently discard.
- For `condition-hide=Average`, a property named exactly `Average`.

The `NA` token, first-column gene IDs, required header, maximum missingness
ceiling of 0.2, and `direction-covar` syntax follow the MAGMA manual linked from
the [official MAGMA documentation](https://cncr.nl/research/magma/).

## Command

```console
postgwas magmacovar --magma-gene-results-file PATH --covariates PATH [options]
```

## Minimal example

Run the default two-sided marginal tests:

```console
postgwas magmacovar \
  --magma-gene-results-file STUDY.genes.raw \
  --covariates gene_covariates.tsv \
  --dataset-id STUDY \
  --output-directory results
```

## Full example

### Original FLAMES

Create the input documented by original FLAMES from its 54-specific-tissue
table. Omitting the two model options retains the configured MAGMA-equivalent
marginal two-sided behavior:

```console
postgwas magmacovar \
  --magma-gene-results-file STUDY.genes.raw \
  --covariates reference/GTEx/gtex_v8_ts_avg_log2TPM.txt \
  --dataset-id STUDY \
  --output-directory results
```

### FUMA-compatible tissue specificity

Run the distinct tissue-specific model used by FUMA's FLAMES wrapper:

```console
postgwas magmacovar \
  --magma-gene-results-file STUDY.genes.raw \
  --covariates reference/GTEx/gtex_v8_ts_general_avg_log2TPM.txt \
  --covariate-model condition-hide=Average \
  --covariate-direction greater \
  --dataset-id STUDY \
  --output-directory results
```

This produces the MAGMA model equivalent to:

```text
--model condition-hide=Average direction-covar=greater
```

FUMA's current FLAMES wrapper consumes the result made from its
30-general-tissue GTEx v8 table, `gtex_v8_ts_general_avg_log2TPM.txt`. This is
not the command printed in the original FLAMES README. The bundled FLAMES
example `.gsa.out` nevertheless records `CONDITIONED_HIDDEN = Average` and a
positive one-sided covariate test, so the bundled result and README command are
internally inconsistent.

## Parameters

`modules.magmacovar` controls:

- `model`: zero or more MAGMA `--model` modifiers; an empty list requests
  separate marginal tests for every property;
- `model_use_cases`: the model, direction, label, and scientific question shown
  in the CLI decision guide;
- `model_placeholders`: definitions of the property, internal-variable, and file
  placeholders shown by the CLI;
- `model_modifiers`: the schema-validated catalogue of modifier syntax,
  published GWAS examples, and example commands shown by the CLI;
- `direction`: `two-sided`, `greater`, or `smaller`;
- `minimum_genes`: the minimum required input overlap and output `NGENES`;
- `input.gene_results_file` and `input.covariates_file`;
- `input.missing_values`: `drop`, `median`, or `mean`;
- `input.maximum_missing_fraction`: between 0 and MAGMA's maximum 0.2;
- `input.missing_genes`: `drop` or `fill`;
- `output_layout`: output, log, resolved-configuration, completion, and staging
  paths.

Multiple valid modifiers can follow `--covariate-model`. They are passed after
one `--model` flag and followed by the explicit
`direction-covar=<direction>` modifier. The schema-validated YAML is the source
of defaults; CLI values only override matching keys. An exported FLAMES pipeline
configuration retains the original README's marginal two-sided defaults; users
can explicitly select the FUMA/bundled-example model instead.

### Reading the advanced model help

`<PROPERTY>` means an exact, case-sensitive header after the leading gene-ID
column in the covariate table; it is not a literal MAGMA value. `<INTERNAL>` is
restricted to `size`/`genesize`, `density`, `mac`, and `sampsize`/`N`. `<PATH>`
points to a whitespace-separated property list, except with `joint`, where each
line defines one multivariable model. Example property names in CLI help are
illustrative and must be replaced by headers present in the selected table.

The research examples translate published analysis designs into current MAGMA
v1.10 syntax; a cited article may describe the regression design without
printing the literal command. They cover four established GWAS uses:

- retaining MAGMA's standard technical correction used in published
  gene-property analyses;
- screening many cell-type properties and restricting conditional follow-up to
  significant results;
- testing FUMA-style tissue or cell-type specificity while adjusting for
  average expression;
- using conditional, residualized, joint, or pairwise models to distinguish
  correlated properties and likely independent associations.

The current gene-property-only interface accepts every documented MAGMA v1.10
modifier that does not require a gene-set annotation: `analyse`, `correct`,
`condition`, `condition-hide`, `condition-residualize`, `joint`, and
`joint-pairs`. Their complete accepted forms and explanations are stored under
`model_modifiers` in the canonical YAML and displayed in
`postgwas magmacovar --help`. Use `--covariate-direction` for MAGMA's
`direction-covar` modifier. Interaction modifiers and gene-set-only modifiers
are deliberately excluded because PostGWAS MAGMAcovar supplies `--gene-covar`
but not `--set-annot`; MAGMA does not support interactions between two
continuous gene properties.

## Processing steps

PostGWAS resolves configuration, validates the executable and both inputs,
checks gene overlap, missingness, and property variation, writes to an isolated
staging prefix, executes MAGMA without shell interpolation, validates the
`.gsa.out` schema and statistics, publishes the result and native log, and writes
the checksummed completion manifest last.

## Outputs

The primary result is `<output>/<dataset>.gsa.out`. MAGMA's native
`<dataset>.log`, the canonical PostGWAS log, resolved configuration, and a
checksummed completion manifest are also retained. Output is published only
after MAGMA exits successfully and the result contains valid `COVAR` rows with
finite statistics and sufficient `NGENES`.

## QC and logs

Review the exact model and direction, MAGMA version and command, gene counts and
overlap, missingness for every property, the number of tested properties,
minimum and maximum result `NGENES`, completion status, and input/output
fingerprints. For tissue-specific analysis, confirm the input table and result
use the intended 30-general- or 54-specific-tissue hypothesis family.

## Interpretation

A property coefficient describes association with gene-level GWAS signal under
the selected MAGMA model and direction. It does not establish mediation or
causality. Review property scaling, missingness, gene overlap, the prespecified
hypothesis family, and multiple-testing correction before interpretation.

## Common problems

Execution stops for a missing MAGMA executable, malformed `.genes.raw` metadata,
duplicate or incompatible gene identifiers, invalid numeric values, insufficient
overlap, excessive missingness, constant properties, or missing output rows.
`--resume` reuses output only when every fingerprint matches; `--overwrite`
replaces only artifacts owned by the configured prefix.

## Limitations

MAGMA text outputs do not encode enough resource metadata to infer the genome
build or prove that gene and covariate identifiers are biologically compatible.
Direct FLAMES mode also cannot recover the MAGMA model from an existing
`.gsa.out`; retain and review its native MAGMA log or equivalent provenance.

## Scientific references

- [Official MAGMA documentation and v1.10 downloads](https://cncr.nl/research/magma/)
- [de Leeuw et al. 2015, MAGMA](https://doi.org/10.1371/journal.pcbi.1004219)
- [de Leeuw et al. 2018, conditional and interaction analysis](https://doi.org/10.1038/s41467-018-06022-6)
- [Watanabe et al. 2019, MAGMA cell-type specificity](https://doi.org/10.1038/s41467-019-11181-1)
- [Duncan et al. 2025, pairwise conditional brain-cell analysis](https://doi.org/10.1038/s41593-024-01834-w)
- [FUMA MAGMA tissue-expression command](https://github.com/vufuma/FUMA-webapp/blob/0c0259b7ed5e6d15369e978530b7ba56b5bd0437/scripts/magma/magma.py#L118-L126)
- [FUMA tissue-expression dataset choices](https://github.com/vufuma/FUMA-webapp/blob/0c0259b7ed5e6d15369e978530b7ba56b5bd0437/resources/views/snp2gene/newjob_components/_magma.blade.php#L61-L87)
- [Original FLAMES MAGMA command](https://github.com/Marijn-Schipper/FLAMES/blob/d411689cfe617a5e8bc719e457dfba05e3f7bd3e/README.md#2-run-magma-tissue-type-analysis-using-your-magma-z-scores-on-the-preformatted-gtex-tissue-expression-file)
- [Bundled FLAMES MAGMA result metadata](https://github.com/Marijn-Schipper/FLAMES/blob/d411689cfe617a5e8bc719e457dfba05e3f7bd3e/example_data/magma_exp_gtex_v8_ts_avg_log2TPM.txt.gsa.out#L1-L4)
