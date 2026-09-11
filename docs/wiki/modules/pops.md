# PoPS

## Purpose

PoPS prioritizes genes by learning genome-wide relationships between MAGMA gene
association scores and a large matrix of gene features. PostGWAS uses the
published upstream PoPS v0.2 algorithm without changing its statistical
implementation.

## What the analysis does

PostGWAS validates and aligns MAGMA gene Z-scores (`ZSTAT`) or custom gene
scores, the PoPS gene
annotation, and every feature-matrix chunk. It then invokes the published
upstream algorithm, verifies its required outputs, and publishes them with
resolved configuration and completion metadata.

## When to use it

Use PoPS after MAGMA gene association when compatible feature matrices are
available and gene identifiers and genome context match across all resources.

## Input requirements

The MAGMA results, PoPS gene annotation, and feature matrices must use compatible
gene identifiers and the same genome build. PostGWAS requires the build to be
declared explicitly because these files do not carry enough metadata to infer it
safely. The service validates gene overlap, duplicate identifiers, finite scores,
feature-chunk dimensions, feature-name uniqueness, chromosome selections, and
finite numeric matrix values, optional feature lists, and optional gene-score
covariance before starting PoPS.

For MAGMA gene Z-scores, preflight also requires `.genes.out` and `.genes.raw` to
contain the same unique genes in the same order, because the raw file supplies
gene covariates and covariance metadata in that order. Every feature-row gene
must occur in the PoPS annotation because it can receive a prediction, and every
gene with an input score must occur in both the annotation and feature rows. Extra
annotation genes are harmless and are allowed. Incompatibility is a hard failure
reported with counts, example identifiers, and source paths under the default
`strict` policy; PostGWAS never silently intersects or discards genes with input
scores.
Shared MAGMA and annotation genes must also have identical chromosome labels;
intersection cannot safely repair a chromosome disagreement.

The MAGMA raw-file preflight also checks the technical-covariate fields used by
PoPS (`NSNPS`, `NPARAM`, and `MAC`) for finite positive values and reconstructs
every chromosome covariance block to verify finite symmetric structure before
the PoPS model starts.

For MAGMA inputs only, `--gene-universe-policy intersect` explicitly restricts
the MAGMA Z-score gene universe to genes present in both the PoPS annotation and
feature rows. This is an opt-in scientific transformation. Original MAGMA files
remain
unchanged. PostGWAS reconstructs covariance-preserving principal submatrices and
publishes separate compatible and excluded `.genes.out`/`.genes.raw` pairs, a
gene-level audit table, and a structured compatibility report. Custom gene-score
files must be aligned by the user and do not support this policy.

PoPS scores rank genes by learned feature similarity to genome-wide association
patterns. They are not calibrated probabilities that a gene is causal.

The completion screen reports the number of genes receiving finite PoPS scores,
the scored genes with and without input gene scores, the genes used to fit the
PoPS model, the genes with input scores not used for model fitting, the selected
feature count, and a configurable top-ranked gene list. The official program
calls this stage “Training” in its source code: it constructs `X_train` and
`Y_train` and estimates regression coefficients. This is model fitting within
the PoPS analysis, not training a separate externally supplied machine-learning
model. PoPS does not define a universal
significance cutoff or a statistically significant gene count: scores are
relative rankings rather than p-values. A top-ranked set should therefore be
prespecified for reporting and supported with independent genetic or functional
evidence. The screen flags incomplete input-score gene or chromosome coverage
because such omissions can make a technically completed model scientifically
unsuitable for genome-wide interpretation.

- A declared `GRCh37` or `GRCh38` genome build.
- A feature prefix with the configured row file and every declared matrix and
  column chunk.
- A gene annotation containing the configured gene-ID, chromosome, and TSS
  columns.
- Either a MAGMA prefix with non-empty `.genes.out` and `.genes.raw` files, or a
  custom gene-score table supplied through `--target-score-file` (the option
  name is retained for compatibility).
- Dataset ID and output directory, supplied by CLI or run configuration.

Pipeline mode supplies the MAGMA prefix from the validated preceding MAGMA step.
It also passes the validated annotated MAGMA gene-result artifact directly to
PoPS for result enrichment; no path is guessed from filenames. PoPS feature
resources are preflighted before any pipeline analysis begins.

PoPS does **not** consume the output of `postgwas magmacovar`. Its
`--use-magma-covariates` option refers to technical gene covariates reconstructed
from MAGMA `.genes.raw` metadata, including gene size, gene density, and inverse
minor-allele-count terms. It does not refer to a tissue-expression `.gsa.out`
file or to `--covariate-model`.

For a standard MAGMA-backed PoPS run, supply:

- `--magma-association-prefix`, shared by `.genes.out` and `.genes.raw`;
- optionally, `--magma-annotated-results-file`, which enriches the integrated
  output with gene symbols, gene intervals, and corrected MAGMA p-values;
- `--feature-matrix-prefix` and the exact `--feature-matrix-chunks` count;
- `--pops-gene-location-file`, using the same Ensembl gene identifiers and
  genome build as MAGMA and the feature rows;
- optionally, `--control-features-file`. The official PoPS example supplies a
  control-feature list, and FUMA's FLAMES integration uses 116 feature chunks
  plus its FUMA-compatible control list.

## Command

```text
postgwas pops --magma-association-prefix PREFIX [options]
postgwas pipeline --modules pops --help
```

## Direct mode

### Existing MAGMA results

```console
postgwas pops \
  --magma-association-prefix magma/STUDY \
  --magma-annotated-results-file magma/STUDY_magma_genes_annotated.tsv \
  --feature-matrix-prefix reference/pops/features_munged/pops_features \
  --feature-matrix-chunks 2 \
  --pops-gene-location-file reference/pops/gene_annot.tsv \
  --control-features-file reference/pops/control.features \
  --genome-build GRCh37 \
  --dataset-id STUDY \
  --output-directory results
```

### Explicit model settings

```console
postgwas pops \
  --magma-association-prefix magma/STUDY \
  --feature-matrix-prefix reference/FLAMES/pops_features_full_FUMA_compatible/features_munged/pops_features \
  --feature-matrix-chunks 116 \
  --pops-gene-location-file reference/FLAMES/pops_features_full_FUMA_compatible/gene_annots.txt \
  --control-features-file reference/FLAMES/pops_features_full_FUMA_compatible/control.features \
  --genome-build GRCh37 \
  --use-magma-covariates \
  --use-magma-error-covariance \
  --feature-selection-p-cutoff 0.05 \
  --remove-hla-during-feature-selection \
  --remove-hla-during-training \
  --method ridge \
  --seed 10 \
  --dataset-id STUDY \
  --output-directory results
```

### Custom gene scores

This direct-only alternative replaces the MAGMA input; do not supply a MAGMA
prefix with it. The packaged schema is a tab-delimited table with `ENSGID` and
`Score` headers, one unique compatible gene per row, and finite numeric scores.
For example, the following illustrates the columns, not a sufficient training
dataset:

```text
ENSGID	Score
ENSG00000139618	2.1
ENSG00000141510	1.7
```

Supply the full prespecified gene-score dataset aligned to the annotation and
feature rows, not a significant-gene-only list. Custom scores cannot use the
MAGMA `intersect` policy. Optional `--target-covariates-file` uses the same
unique `ENSGID` set plus finite numeric covariates; optional
`--target-error-covariance-file` is a finite symmetric `.npy`/sparse `.npz`
matrix in the score table's gene order. Their derivation and score scale need
independent justification; an arbitrary score is not equivalent to MAGMA ZSTAT.

```console
postgwas pops \
  --target-score-file target_scores.tsv \
  --feature-matrix-prefix reference/pops/features_munged/pops_features \
  --feature-matrix-chunks 2 \
  --pops-gene-location-file reference/pops/gene_annot.tsv \
  --genome-build GRCh37 \
  --method ridge \
  --dataset-id STUDY \
  --output-directory results/pops_custom
```

### Alternative regression models

`ridge` is the packaged default; `lasso` and `linreg` are also supported.
Prespecify the model rather than choosing whichever yields preferred genes.
These alternatives use the same input validation, feature-selection and
chromosome policies; changing the fit can change rankings.

```console
postgwas pops \
  --magma-association-prefix magma/STUDY \
  --feature-matrix-prefix reference/pops/features_munged/pops_features \
  --feature-matrix-chunks 2 \
  --pops-gene-location-file reference/pops/gene_annot.tsv \
  --genome-build GRCh37 \
  --method lasso \
  --dataset-id STUDY \
  --output-directory results/pops_lasso
```

```console
postgwas pops \
  --magma-association-prefix magma/STUDY \
  --feature-matrix-prefix reference/pops/features_munged/pops_features \
  --feature-matrix-chunks 2 \
  --pops-gene-location-file reference/pops/gene_annot.tsv \
  --genome-build GRCh37 \
  --method linreg \
  --dataset-id STUDY \
  --output-directory results/pops_linreg
```

### Optional MAGMA gene-universe intersection

Use this only after reviewing why the MAGMA and PoPS gene releases differ:

```console
postgwas pops \
  --magma-association-prefix magma/STUDY \
  --feature-matrix-prefix reference/pops/features_munged/pops_features \
  --feature-matrix-chunks 116 \
  --pops-gene-location-file reference/pops/gene_annot.tsv \
  --gene-universe-policy intersect \
  --genome-build GRCh37 \
  --dataset-id STUDY \
  --output-directory results
```

## Pipeline mode

Use a strand-aware MAGMA gene-location file from the same Ensembl gene universe
as the PoPS annotation and feature rows. An NCBI/Entrez location file is not an
alternative for an Ensembl feature bundle. The following GRCh37/EUR example
assumes a 116-chunk feature bundle. Replace the file paths, chunk count, and the
uppercase source-metadata placeholders with the exact values from your resource
manifest; the source URL must be HTTPS. Metadata declarations do not remap IDs.

```console
postgwas pipeline \
  --modules pops \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --magma-ld-reference reference/1000G_EUR \
  --gene-location-file reference/pops/PoPS_GRCh37_strand_aware.loc \
  --magma-positional-gene-id-type ensembl \
  --magma-positional-source-name GENE_LOCATION_SOURCE \
  --magma-positional-source-version SOURCE_RELEASE \
  --magma-positional-source-url https://example.org/GENE_LOCATION_SOURCE_RECORD \
  --magma-positional-context GENE_LOCATION_CONTEXT \
  --feature-matrix-prefix reference/pops/features_munged/pops_features \
  --feature-matrix-chunks 116 \
  --pops-gene-location-file reference/pops/GRCh37_gene_annot.tsv \
  --pops-genome-build GRCh37 \
  --dataset-id STUDY \
  --output-directory results
```

The pipeline supplies MAGMA's validated `.genes.out`/`.genes.raw` pair; do not
provide `--magma-association-prefix` in pipeline mode. The direct build option
is `--genome-build`; the PoPS-specific pipeline option is `--pops-genome-build`.
For non-default MAGMA builds or mappings, also provide compatible MAGMA settings
as described in the [MAGMA guide](magma.md#gene-identifiers-and-source-declarations).

The preceding pipeline uses the default `ridge` model. Complete alternatives
for the same GRCh37/EUR, 116-chunk resource contract follow. Replace all source
metadata and resource placeholders exactly as above.

### lasso regression

```console
postgwas pipeline \
  --modules pops \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --magma-ld-reference reference/1000G_EUR \
  --gene-location-file reference/pops/PoPS_GRCh37_strand_aware.loc \
  --magma-positional-gene-id-type ensembl \
  --magma-positional-source-name GENE_LOCATION_SOURCE \
  --magma-positional-source-version SOURCE_RELEASE \
  --magma-positional-source-url https://example.org/GENE_LOCATION_SOURCE_RECORD \
  --magma-positional-context GENE_LOCATION_CONTEXT \
  --feature-matrix-prefix reference/pops/features_munged/pops_features \
  --feature-matrix-chunks 116 \
  --pops-gene-location-file reference/pops/GRCh37_gene_annot.tsv \
  --pops-genome-build GRCh37 \
  --method lasso \
  --dataset-id STUDY \
  --output-directory results/pops_lasso
```

### linreg regression

```console
postgwas pipeline \
  --modules pops \
  --vcf STUDY_GRCh37_merged.vcf.gz \
  --magma-ld-reference reference/1000G_EUR \
  --gene-location-file reference/pops/PoPS_GRCh37_strand_aware.loc \
  --magma-positional-gene-id-type ensembl \
  --magma-positional-source-name GENE_LOCATION_SOURCE \
  --magma-positional-source-version SOURCE_RELEASE \
  --magma-positional-source-url https://example.org/GENE_LOCATION_SOURCE_RECORD \
  --magma-positional-context GENE_LOCATION_CONTEXT \
  --feature-matrix-prefix reference/pops/features_munged/pops_features \
  --feature-matrix-chunks 116 \
  --pops-gene-location-file reference/pops/GRCh37_gene_annot.tsv \
  --pops-genome-build GRCh37 \
  --method linreg \
  --dataset-id STUDY \
  --output-directory results/pops_linreg
```

## Parameters

Export the canonical YAML and edit it for reusable runs:

```console
postgwas config export \
  --module pops \
  --style full \
  --output pops.yaml

postgwas pops --run-config pops.yaml
```

Explicit CLI values override matching YAML values. All omitted values come from
the schema-validated YAML; the CLI does not maintain a second set of defaults.
See the generated [Configuration Defaults](../reference/configuration-defaults.md)
for every input schema, model, feature-selection, output, and continuation key.

## Processing steps

The direct PoPS screen and canonical log report these stages in order:

1. validate the MAGMA gene Z-scores or custom gene scores and optional
   gene-score covariates and covariance;
2. validate the PoPS gene annotation, including its gene-ID, chromosome, and
   transcription-start-site (TSS) columns;
3. validate the feature-row file and every configured column/matrix chunk,
   including dimensions, numeric data type, finite values, and globally unique
   feature names;
4. validate optional feature-subset and control-feature lists against the full
   matrix feature universe, and ensure a supplied subset does not silently
   remove requested controls;
5. compare input-score, annotation, and feature-row gene identifiers and apply
   the configured gene-universe policy;
6. load the gene scores used to fit PoPS—MAGMA gene Z-scores (`ZSTAT`) for a
   MAGMA-backed run or the configured score column for a custom gene-score run;
7. adjust those PoPS fitting scores for configured covariates;
8. test and select predictive features;
9. fit the prediction model;
10. calculate genome-wide PoPS scores; and
11. validate the official PoPS outputs, build the integrated full-union TSV and
    HTML report, and publish all required results atomically.

A PoPS-only pipeline first displays the eight gene-only MAGMA prerequisite
stages—beginning with validation of the input summary-statistics VCF, PLINK LD
reference, and MAGMA gene-location reference—then the eleven PoPS stages above.
Pathway stages are not included because PoPS consumes MAGMA gene-association
results, not MAGMA pathway results.

PoPS uses the configured TSS column to apply its chromosome 6 HLA-region
exclusions during covariate projection, feature selection, and model training
when those policies are enabled. The annotation is not used to redefine MAGMA
gene intervals.

PostGWAS writes the resolved configuration and runs upstream PoPS in an isolated
staging directory. Required output files are checked for existence and content
before publication. A completion manifest is written last; outputs without this
manifest are incomplete and are never resumed.

## Outputs

With the default output prefix `<output>/<dataset>_pops`, the published upstream
files are `.preds`, `.coefs`, `.marginals`, and `.log`. Optional `.traindata` and
`.matdata` files are published when matrix saving is enabled. PostGWAS also
writes a structured service log, resolved configuration, and completion manifest.
The prediction table contains `ENSGID`, `PoPS_Score`, the input gene score `Y`,
and the upstream feature-selection and model-fitting flags. Every feature-matrix
row receives a PoPS score, including a feature-row gene whose `Y` is missing.
Such a gene was not used to select features or fit coefficients, but it may
still have a high or low PoPS score because the fitted coefficients are applied
to all feature rows. PostGWAS preserves this official behavior and identifies
those genes explicitly instead of silently removing them.

Two additional outputs are always published:

- `<dataset>_pops.integrated_gene_results.tsv`: a full outer union of genes in
  the official `.preds` file and the original MAGMA or custom gene-score table;
- `<dataset>_pops.integrated_gene_results.html`: a standalone searchable,
  sortable, paginated view of the same records, with a TSV download link and
  category counts.

The integrated table does not replace or modify any authoritative source file.
Its leading columns are arranged for interpretation: gene identity and PoPS
coordinates, PoPS score and rank, analysis status, and the principal MAGMA
association statistics. Detailed input-score alignment, model-use, availability,
and reference-coordinate fields follow as an audit trail. It records source
availability, the original input row number and gene score, PoPS score and rank,
the `Y` and optional covariate-adjusted `Y_proj` values written
by PoPS, covariate-projection/feature-selection/model-fitting flags, compatibility
decisions, PoPS annotation chromosome/TSS, and available MAGMA statistics and
annotations. When an annotated MAGMA table is supplied, its gene set,
chromosomes, and Z statistics must agree with `.genes.out`/`.genes.raw` before
PoPS starts. Without it, fields available directly from `.genes.out` are still
included and enrichment-only fields are `NA`; custom gene-score runs leave all
MAGMA-specific fields `NA`.

`gene_analysis_status` gives every union row one unambiguous category:

- `scored_and_used_for_model_fitting`;
- `scored_with_target_not_used_for_model_fitting` (for example, a gene removed
  from fitting by the configured chromosome or HLA policy);
- `scored_without_input_target_score`;
- `input_target_excluded_from_pops` (possible only with explicit MAGMA
  intersection).

The distinction between `input_target_score` and `pops_target_score` is
intentional. These machine-readable column names are retained for compatibility.
The first comes from the original MAGMA/custom gene-score table; the second is
the `Y` value written into `.preds` by upstream PoPS after its validated gene
alignment. PostGWAS requires the two to agree for every retained gene with an
input score.

With `gene_universe_policy: intersect`, six additional owned outputs are
published:

- `.compatible.genes.out` and `.compatible.genes.raw`: the aligned MAGMA pair
  used by PoPS;
- `.excluded.genes.out` and `.excluded.genes.raw`: excluded genes and their
  covariance-preserving principal submatrix;
- `.gene_compatibility.tsv`: one decision row per original MAGMA gene;
- `.gene_compatibility.yaml`: total and per-chromosome counts, percentages,
  exclusion reasons, input/output paths, covariance handling, and the scientific
  effect of filtering.

NumPy `.npy` and SciPy sparse `.npz` custom covariance inputs are validated.
NumPy covariance is converted in the staging directory to the sparse container
expected by upstream PoPS; scientific values and ordering are unchanged.

## QC and logs

Review the declared genome build, input gene-score structure, gene-annotation
columns and TSS range, shared-gene counts, optional control-list compatibility,
chromosome/HLA policies, method and seed, upstream log, structured service log,
resolved configuration, and completion manifest. PostGWAS validates every
feature chunk's dimensions, numeric type, and finite values. To keep a 116-chunk
run readable, the terminal shows only the configured number of successful file
sections and counts the remainder; the canonical input-validation YAML retains
every column and matrix file with its individual checks. For intersection runs,
also review
retained/excluded counts and percentages, missing-resource categories, every
gene-level decision, and per-chromosome statistics.

## Interpretation

PoPS scores are relative gene-prioritization scores learned from polygenic gene
features. They are not causal probabilities and should be interpreted alongside
locus-based and other orthogonal evidence.

## Common problems

Invalid or incompatible resources fail before model fitting with an actionable
message. Failed upstream runs remain isolated below the configured staging
directory and are not presented as completed outputs. With resume enabled, only
a complete output set carrying its completion manifest is reused. Use
`--overwrite` to replace an existing result intentionally.

## Limitations

PostGWAS cannot infer the genome build or ancestry from legacy PoPS matrices or
MAGMA text outputs, so compatibility remains a user-declared scientific
responsibility. The unchanged upstream implementation constructs dense
covariance matrices and can require substantial memory for large gene sets.

## Scientific references

- [Weeks et al. 2023, PoPS](https://doi.org/10.1038/s41588-023-01443-6)
- [Official PoPS v0.2 implementation](https://github.com/FinucaneLab/pops)
- [FUMA FLAMES PoPS invocation](https://github.com/vufuma/FUMA-webapp/blob/0c0259b7ed5e6d15369e978530b7ba56b5bd0437/scripts/flames/run_flames.py#L198-L221)
