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

PostGWAS applies the configured multiple-testing methods globally across all
validated `COVAR` rows from one MAGMA run. The packaged methods are Bonferroni,
Šidák, Holm, and Benjamini-Hochberg FDR, with Bonferroni declared as the primary
method at an inclusive adjusted-p threshold of 0.05. This primary choice matches
the Bonferroni correction used for the 461-cell-property screen in Duncan et al.
(2025). These settings are YAML-controlled and may be changed only as a
prespecified analysis decision.

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

### Pre-MAGMA missingness guard

MAGMA removes a gene property when its missingness is greater than `max-miss`.
PostGWAS calculates the missing fraction for every property before starting
MAGMA and raises an error that names any failing property. This prevents a
silently removed tissue or other property from being mistaken for a test that
completed without a result. A property exactly equal to `max-miss` is accepted;
only a value strictly greater than the threshold fails.

The denominator follows MAGMA's two missing-gene policies. With
`missing-genes=fill`, it is the number of eligible genes in `.genes.raw`, and
genes absent from the covariate table are added to the missing count. Otherwise,
it is the number of gene IDs overlapping the two files, because absent genes are
dropped. Reversing these denominators could accept a property MAGMA would remove
or reject a property MAGMA would retain.

For the direction, PostGWAS emits the property-specific
`direction-covar=<DIRECTION>` modifier rather than the blanket `direction=` or
gene-set-specific `direction-sets=` modifier. The accepted values use MAGMA's
vocabulary: `two-sided`, `greater`, and `smaller` (`less` is not accepted).

The three missing-data settings can be supplied through YAML or overridden on
either the direct or pipeline command line:

| Purpose | CLI option | YAML key | Accepted values | Default |
| --- | --- | --- | --- | --- |
| Handle `NA` values after the property passes `max-miss` | `--covariate-missing-values ACTION` | `modules.magmacovar.input.missing_values` | `drop`, `median`, `mean` | `median` |
| Set the per-property missingness threshold | `--covariate-max-miss FRACTION` | `modules.magmacovar.input.maximum_missing_fraction` | `0` through `0.2` | `0.05` |
| Handle eligible genes absent from the covariate table | `--covariate-missing-genes ACTION` | `modules.magmacovar.input.missing_genes` | `drop`, `fill` | `drop` |

Correcting or removing a sparse property is preferred. Increase
`--covariate-max-miss` only with a documented data-quality justification, and
never above MAGMA's limit of `0.2`. Changing the missing-value replacement
method does not bypass the missingness threshold.

## Command

```text
postgwas magmacovar --magma-gene-results-file PATH --covariates PATH [options]
```

## Direct mode

### Marginal gene-property tests

Run the default two-sided marginal tests:

```console
postgwas magmacovar \
  --magma-gene-results-file STUDY.genes.raw \
  --covariates gene_covariates.tsv \
  --dataset-id STUDY \
  --output-directory results
```

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

The following are distinct model designs. Replace property names with exact
headers in your covariate table and prespecify the tested family. For the
file-based joint example, `property_models.txt` contains one whitespace-separated
list of properties per model, one model per line. It is not a list of genes.
`analyse` restricts the tested properties; it is not a significance filter.
Technical corrections remain at the configured MAGMA default unless explicitly
changed with `correct`.

### Conditional properties

```console
postgwas magmacovar \
  --magma-gene-results-file STUDY.genes.raw \
  --covariates gene_covariates.tsv \
  --covariate-model analyse=list,CellType_A condition=CellType_B \
  --dataset-id STUDY \
  --output-directory results
```

### Residualized outcome

```console
postgwas magmacovar \
  --magma-gene-results-file STUDY.genes.raw \
  --covariates gene_covariates.tsv \
  --covariate-model analyse=list,BiologicalScore condition-residualize=TechnicalScore \
  --dataset-id STUDY \
  --output-directory results
```

### Joint models from a file

```console
postgwas magmacovar \
  --magma-gene-results-file STUDY.genes.raw \
  --covariates gene_covariates.tsv \
  --covariate-model joint=property_models.txt \
  --dataset-id STUDY \
  --output-directory results
```

### Pairwise joint models

```console
postgwas magmacovar \
  --magma-gene-results-file STUDY.genes.raw \
  --covariates gene_covariates.tsv \
  --covariate-model analyse=list,CellType_A,CellType_B joint-pairs \
  --dataset-id STUDY \
  --output-directory results
```

## Pipeline mode

### Marginal gene-property tests

When `magmacovar` is the requested pipeline target, PostGWAS starts from the
summary-statistics GWAS-VCF and runs its MAGMA dependency before gene-property
analysis. The pipeline command therefore accepts the GWAS-VCF, PLINK LD
reference prefix, gene-location reference, and gene-property covariate table:

```console
postgwas pipeline \
  --modules magmacovar \
  --vcf study.vcf.gz \
  --magma-ld-reference reference/1000G_EUR \
  --gene-location-file reference/NCBI37.3.gene.loc \
  --covariates gene_covariates.tsv \
  --dataset-id STUDY \
  --output-directory results
```

The LD, gene-location, build, population, and mapping settings remain owned by
the schema-validated `modules.magma` configuration. MAGMAcovar reuses that
resolved configuration and the validated `.genes.raw` handoff; it does not
define a second copy of those reference settings. The MAGMA dependency is
gene-only for this target: it does not validate a pathway file, align pathway
identifiers, run competitive pathway association, or publish pathway results.
Use the standalone `magma` pipeline target when pathway association is wanted.
Explicit pathway controls, including `--gene-set-file` and its identifier-
compatibility settings, are rejected on a MAGMAcovar-only command so a
requested pathway analysis is never silently discarded. The pipeline help is
kept concise; the complete advanced gene-property model reference remains in
`postgwas magmacovar --help`.

For the canonical single positional-MAGMA mapping (the packaged default), the
shared screen progress and durable screen log show one ordered 12-stage account:

1. Stages 1–3 validate the GWAS-VCF, PLINK BED/BIM/FAM reference, BIM identifier
   convention, and gene-location file.
2. Stages 4–8 prepare the variant inputs, create the SNP-to-gene annotation,
   run and validate gene association, and publish the validated MAGMA outputs.
3. Stages 9–12 validate the native `.genes.raw` handoff, validate the covariate
   table and its gene overlap and missingness against the MAGMA genes, run the
   configured gene-property model, and validate and correct the `COVAR` results.
   Output publication and checkpoint creation remain internal completion checks,
   so they do not add a user-facing progress stage or outcome block.

Each completed stage displays its observed counts, decisions, filenames, and
validation status. The same evidence is retained in the screen transcript; the
MAGMA canonical log records the upstream reference and association stages, and
the MAGMAcovar canonical log records each of its four stages and outcomes. If a
pathway file is configured elsewhere in reusable MAGMA settings, it is outside
the MAGMAcovar gene-only execution contract and is not validated or used in a
MAGMA command.

These examples use the packaged GRCh37/EUR positional MAGMA reference
contract and an Entrez-compatible covariate table. For a different primary gene
identifier system, supply a matching MAGMA reference and its
[source declarations](magma.md#gene-identifiers-and-source-declarations).
The model meanings and property-file requirements are the same as in direct mode.

### Tissue specificity

```console
postgwas pipeline \
  --modules magmacovar \
  --vcf study_GRCh37.vcf.gz \
  --magma-ld-reference reference/1000G_EUR \
  --gene-location-file reference/NCBI37.3.gene.loc \
  --covariates gene_covariates.tsv \
  --covariate-model condition-hide=Average \
  --covariate-direction greater \
  --dataset-id STUDY \
  --output-directory results
```

### Conditional properties

```console
postgwas pipeline \
  --modules magmacovar \
  --vcf study_GRCh37.vcf.gz \
  --magma-ld-reference reference/1000G_EUR \
  --gene-location-file reference/NCBI37.3.gene.loc \
  --covariates gene_covariates.tsv \
  --covariate-model analyse=list,CellType_A condition=CellType_B \
  --covariate-direction two-sided \
  --dataset-id STUDY \
  --output-directory results
```

### Residualized outcome

```console
postgwas pipeline \
  --modules magmacovar \
  --vcf study_GRCh37.vcf.gz \
  --magma-ld-reference reference/1000G_EUR \
  --gene-location-file reference/NCBI37.3.gene.loc \
  --covariates gene_covariates.tsv \
  --covariate-model analyse=list,BiologicalScore condition-residualize=TechnicalScore \
  --covariate-direction two-sided \
  --dataset-id STUDY \
  --output-directory results
```

### Joint models from a file

```console
postgwas pipeline \
  --modules magmacovar \
  --vcf study_GRCh37.vcf.gz \
  --magma-ld-reference reference/1000G_EUR \
  --gene-location-file reference/NCBI37.3.gene.loc \
  --covariates gene_covariates.tsv \
  --covariate-model joint=property_models.txt \
  --covariate-direction two-sided \
  --dataset-id STUDY \
  --output-directory results
```

### Pairwise joint models

```console
postgwas pipeline \
  --modules magmacovar \
  --vcf study_GRCh37.vcf.gz \
  --magma-ld-reference reference/1000G_EUR \
  --gene-location-file reference/NCBI37.3.gene.loc \
  --covariates gene_covariates.tsv \
  --covariate-model analyse=list,CellType_A,CellType_B joint-pairs \
  --covariate-direction two-sided \
  --dataset-id STUDY \
  --output-directory results
```

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
- `multiple_testing.methods`: ordered global corrections applied across every
  valid `COVAR` result in the run;
- `multiple_testing.primary_method`: method used for the primary significance
  declaration;
- `multiple_testing.significance_threshold`: inclusive adjusted-p threshold for
  that primary declaration;
- `multiple_testing.reporting_method_labels`: terminal labels for the methods;
- `result_schema`: source-statistic, adjusted-p, primary-decision, delimiter,
  and null-value fields for the corrected report;
- `reporting.top_property_count`: number of leading properties shown in the
  main findings after ranking by raw MAGMA `P` (default `5`);
- `reporting.highlight_method`: secondary adjusted-p method highlighted in the
  findings (default `fdr_bh`, displayed as BH-FDR);
- `reporting.p_value_significant_digits`: display precision for raw and
  highlighted adjusted p-values in those findings;
- `reporting.effect_significant_digits`: display precision for standardized
  property coefficients in those findings;
- `output_layout`: output, log, resolved-configuration, completion, and staging
  paths.

Multiple valid modifiers can follow `--covariate-model`. They are passed after
one `--model` flag and followed by the explicit
`direction-covar=<direction>` modifier. The schema-validated YAML is the source
of defaults; CLI values only override matching keys. An exported FLAMES pipeline
configuration retains the original README's marginal two-sided defaults; users
can explicitly select the FUMA/bundled-example model instead.

The corresponding CLI overrides are `--covariate-missing-values`,
`--covariate-max-miss`, and `--covariate-missing-genes`. They affect the same
schema-validated keys listed above and therefore do not introduce separate CLI
defaults.

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
`.gsa.out` schema and statistics, computes each configured global correction,
validates the corrected report against the native rows and configured methods,
publishes both result files and the native log, and writes the checksummed
completion manifest last.

## Outputs

Two scientific result files are published:

- `<output>/<dataset>.gsa.out` is MAGMA's unmodified native result. Its `P`
  column remains raw. This is a MAGMAcovar result, not an extra FLAMES-only
  artifact: PostGWAS retains it as the source for correction validation,
  scientific provenance, and checkpoint reuse. If FLAMES is also requested,
  this native file—not the corrected report—is its downstream handoff because
  FLAMES expects MAGMA's native `VARIABLE` and `P` contract.
- `<output>/<dataset>_magmacovar_corrected.tsv` is the PostGWAS interpretation
  report. It retains `VARIABLE`, `TYPE`, `NGENES`, `BETA`, `BETA_STD`, `SE`, and
  raw `P`, then adds `P_bonferroni_corr`, `P_sidak_corr`, `P_holm_corr`,
  `P_fdr_bh_corr`, `primary_correction_method`, `primary_adjusted_p`, and
  `primary_significant` under the packaged configuration.

MAGMA's native `<dataset>.log`, the canonical PostGWAS log, resolved
configuration, and a checksummed completion manifest are also retained. All
three native/corrected/log artifacts participate in resume identity and output
fingerprinting. Nothing is published until MAGMA exits successfully, the raw
result contains valid `COVAR` rows with finite statistics and sufficient
`NGENES`, and the corrected report exactly matches those rows and the resolved
configuration.

The terminal and screen log summarize the correction family, Bonferroni-primary
and BH-FDR significant-property counts, gene-count range, model, direction, and
the configured number of top properties. Top properties are ordered by the
lowest native MAGMA `P` with property name, standardized coefficient, raw `P`,
and BH-FDR-adjusted p-value. When MAGMA supplies `FULL_NAME`, PostGWAS uses it
instead of the shortened display form in `VARIABLE`, so long tissue or cell-type
names remain complete in both the corrected TSV and the findings. The
`FLAMES handoff` row is displayed only when FLAMES is actually part of the
requested pipeline; it is absent from standalone and MAGMAcovar-only runs.

## QC and logs

Review the exact model and direction, MAGMA version and command, gene counts and
overlap, missingness for every property, the number of tested properties,
minimum and maximum result `NGENES`, the complete correction family, declared
primary method, primary and BH-FDR significant-property counts, top-ranked
properties, completion status, and input/output fingerprints. For
tissue-specific analysis, confirm the input table and result use the intended
30-general- or 54-specific-tissue hypothesis family.

## Interpretation

A property coefficient describes association with gene-level GWAS signal under
the selected MAGMA model and direction. It does not establish mediation or
causality. Review property scaling, missingness, gene overlap, the prespecified
hypothesis family, and multiple-testing correction before interpretation. Use
`primary_adjusted_p` and `primary_significant` for the prespecified primary
claim; the other adjusted columns are secondary sensitivity summaries, not a
menu from which to select the most favorable result after seeing the data.

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
