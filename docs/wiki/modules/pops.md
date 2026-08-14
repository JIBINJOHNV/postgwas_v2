# PoPS

## Purpose

PoPS prioritizes genes by learning genome-wide relationships between MAGMA gene
association scores and a large matrix of gene features. PostGWAS uses the
published upstream PoPS v0.2 algorithm without changing its statistical
implementation.

## What the analysis does

PostGWAS validates and aligns MAGMA or custom target scores, the PoPS gene
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
optional target covariance before starting PoPS.

For MAGMA targets, preflight also requires `.genes.out` and `.genes.raw` to
contain the same unique genes in the same order, because the raw file supplies
gene covariates and covariance metadata in that order. Every feature-row gene
must occur in the PoPS annotation because it can receive a prediction, and every
target-scored gene must occur in both the annotation and feature rows. Extra
annotation genes are harmless and are allowed. Incompatibility is a hard failure
reported with counts, example identifiers, and source paths under the default
`strict` policy; PostGWAS never silently intersects or discards target genes.
Shared MAGMA and annotation genes must also have identical chromosome labels;
intersection cannot safely repair a chromosome disagreement.

For MAGMA inputs only, `--gene-universe-policy intersect` explicitly restricts
the target universe to genes present in both the PoPS annotation and feature
rows. This is an opt-in scientific transformation. Original MAGMA files remain
unchanged. PostGWAS reconstructs covariance-preserving principal submatrices and
publishes separate compatible and excluded `.genes.out`/`.genes.raw` pairs, a
gene-level audit table, and a structured compatibility report. Custom target
files must be aligned by the user and do not support this policy.

PoPS scores rank genes by learned feature similarity to genome-wide association
patterns. They are not calibrated probabilities that a gene is causal.

The completion screen reports the number of genes scored, the number with
MAGMA or custom target scores, the genes used for training, selected feature
count, and a configurable top-ranked gene list. PoPS does not define a universal
significance cutoff or a statistically significant gene count: scores are
relative rankings rather than p-values. A top-ranked set should therefore be
prespecified for reporting and supported with independent genetic or functional
evidence. The screen flags incomplete target-gene or chromosome coverage because
such omissions can make a technically completed model scientifically unsuitable
for genome-wide interpretation.

- A declared `GRCh37` or `GRCh38` genome build.
- A feature prefix with the configured row file and every declared matrix and
  column chunk.
- A gene annotation containing the configured gene-ID, chromosome, and TSS
  columns.
- Either a MAGMA prefix with non-empty `.genes.out` and `.genes.raw` files, or a
  custom target-score table.
- Dataset ID and output directory, supplied by CLI or run configuration.

Pipeline mode supplies the MAGMA prefix from the validated preceding MAGMA step.
PoPS feature resources are preflighted before any pipeline analysis begins.

PoPS does **not** consume the output of `postgwas magmacovar`. Its
`--use-magma-covariates` option refers to technical gene covariates reconstructed
from MAGMA `.genes.raw` metadata, including gene size, gene density, and inverse
minor-allele-count terms. It does not refer to a tissue-expression `.gsa.out`
file or to `--covariate-model`.

For a standard MAGMA-backed PoPS run, supply:

- `--magma-association-prefix`, shared by `.genes.out` and `.genes.raw`;
- `--feature-matrix-prefix` and the exact `--feature-matrix-chunks` count;
- `--pops-gene-location-file`, using the same Ensembl gene identifiers and
  genome build as MAGMA and the feature rows;
- optionally, `--control-features-file`. The official PoPS example supplies a
  control-feature list, and FUMA's FLAMES integration uses 116 feature chunks
  plus its FUMA-compatible control list.

## Command

```console
postgwas pops --magma-association-prefix PREFIX [options]
```

## Minimal example

```console
postgwas pops \
  --magma-association-prefix magma/STUDY \
  --feature-matrix-prefix reference/pops/features_munged/pops_features \
  --feature-matrix-chunks 2 \
  --pops-gene-location-file reference/pops/gene_annot.tsv \
  --control-features-file reference/pops/control.features \
  --genome-build GRCh37 \
  --dataset-id STUDY \
  --output-directory results
```

## Full example

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

PostGWAS validates all resources, writes the resolved configuration, and runs
upstream PoPS in an isolated staging directory. Required output files are checked
for existence and content before publication. A completion manifest is written
last; outputs without this manifest are incomplete and are never resumed.

## Outputs

With the default output prefix `<output>/<dataset>_pops`, the published upstream
files are `.preds`, `.coefs`, `.marginals`, and `.log`. Optional `.traindata` and
`.matdata` files are published when matrix saving is enabled. PostGWAS also
writes a structured service log, resolved configuration, and completion manifest.
The prediction table contains `ENSGID` and `PoPS_Score`.

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

Review the declared genome build, shared-gene counts, feature chunk dimensions,
feature count, chromosome/HLA policies, method and seed, upstream log, structured
service log, resolved configuration, and completion manifest. For intersection
runs, also review retained/excluded counts and percentages, missing-resource
categories, every gene-level decision, and per-chromosome statistics.

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
