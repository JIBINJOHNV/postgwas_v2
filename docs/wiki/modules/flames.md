# FLAMES

## Purpose

FLAMES prioritizes candidate effector genes at fine-mapped GWAS loci by combining
variant-to-gene annotations, MAGMA gene association and tissue evidence, and
PoPS scores. PostGWAS validates and stages the inputs, runs the unchanged
upstream FLAMES program, validates its outputs, and publishes provenance
metadata.

## What the analysis does

PostGWAS consumes the explicit `indexfile.txt` created by its fine-mapping
handoff. It validates every indexed credible set, MAGMA table, PoPS table,
annotation bundle, model file, feature manifest, and runtime dependency before
annotation. Annotation and scoring use the same generated index and execute as
argument arrays without shell interpolation.

The bundled published classifier is calibrated for genes within 750 kb of a
locus. That distance is therefore a fixed model invariant, not an adjustable
analysis default. PostGWAS also requires each input credible set to contain at
least 0.95 cumulative posterior inclusion probability, within the configured
numerical tolerance, and rejects total probability above one.

## When to use it

Use FLAMES after compatible fine-mapping, MAGMA gene association, MAGMA
gene-property analysis, and PoPS scoring have completed. The credible-set
coordinates and all declared module genome builds must agree. MAGMA and PoPS
must use compatible Ensembl gene identifiers; PostGWAS matches identifiers
exactly and does not strip versions or silently remap genes.

## Input requirements

- A PostGWAS fine-mapping interchange directory containing the configured
  `indexfile.txt` and every credible-set file named by that index.
- Credible-set variants matching the configured chromosome-position-allele
  identifier contract, with unique variants and finite PIP values in `[0, 1]`.
- MAGMA `.genes.out`, MAGMA gene-property `.gsa.out`, and PoPS `.preds` files
  containing the configured columns and finite scientific values.
- A prespecified MAGMA gene-property model. Pipeline defaults (`model: []`,
  `direction: two-sided`) reproduce the literal original FLAMES README command.
  FUMA and the bundled FLAMES example instead use `condition-hide=Average` and
  `greater`; PostGWAS accepts that explicitly selected alternative.
- The FLAMES annotation bundle identified in YAML and the packaged or explicitly
  configured model and feature manifest.
- A Python runtime containing every dependency listed in canonical YAML.
- For local VEP, an executable and build-compatible cache. For local CADD, a
  build-compatible bgzip/tabix file, its `.tbi` index, and tabix executable.

## Command

```console
postgwas flames [--run-config PATH] [overrides]
postgwas pipeline --modules flames [pipeline options]
```

## Minimal example

Export and edit the canonical schema-validated configuration:

```console
postgwas config export \
  --module flames \
  --style full \
  --output flames.yaml

postgwas flames --run-config flames.yaml
```

At minimum, set `run.dataset_id`, `run.output_directory`, and these FLAMES
module keys in the exported YAML: `credible_sets_directory`,
`magma_gene_results_file`, `magma_covariate_results_file`, `pops_scores_file`,
and `annotation_resource_directory`.

## Full example

### Direct mode

Run FLAMES from validated existing upstream outputs. Values not shown remain
resolved from canonical YAML. Direct mode validates the `.gsa.out` structure but
cannot recover the MAGMA model from that result file; verify its paired MAGMA log
or other provenance before using it:

```console
postgwas flames \
  --credible-sets-directory finemap/downstream_inputs/flames \
  --magma-gene-results-file magma/STUDY.genes.out \
  --magma-covariate-results-file magma/STUDY.gsa.out \
  --pops-scores-file pops/STUDY.preds \
  --flames-annotation-directory reference/FLAMES/Annotation_data \
  --flames-genome-build GRCh37 \
  --dataset-id STUDY \
  --output-directory results
```

### Pipeline mode

The pipeline creates the fine-mapping, MAGMA, gene-property, and PoPS inputs
consumed by FLAMES:

```console
postgwas pipeline \
  --modules flames \
  --vcf study_GRCh37.vcf.gz \
  --genome-build GRCh37 \
  --ld-region-dir reference/ld_blocks \
  --ld-block-populations EUR \
  --ld-folder reference/pairwise_ld \
  --population EUR \
  --magma-ld-reference reference/1000G_EUR \
  --gene-location-file reference/FUMA/ENSGv102.coding.genes.txt \
  --covariates reference/GTEx/gtex_v8_ts_avg_log2TPM.txt \
  --feature-matrix-prefix reference/FLAMES/pops_features_full_FUMA_compatible/features_munged/pops_features \
  --feature-matrix-chunks 116 \
  --pops-gene-location-file reference/FLAMES/pops_features_full_FUMA_compatible/gene_annots.txt \
  --control-features-file reference/FLAMES/pops_features_full_FUMA_compatible/control.features \
  --pops-genome-build GRCh37 \
  --finemap-method susie \
  --finemap-ld-reference reference/1000G_EUR \
  --flames-annotation-directory reference/FLAMES/Annotation_data \
  --flames-genome-build GRCh37 \
  --plink plink \
  --bcftools bcftools \
  --dataset-id STUDY \
  --output-directory results
```

The FLAMES pipeline defaults to the original README's marginal, two-sided
MAGMACOVAR configuration. `postgwas config export --pipeline flames` records
those values in the generated YAML. No model flags are needed in an all-CLI run
because those values are canonical defaults.

FUMA's current FLAMES wrapper passes the gene-property result for its
30-general-tissue GTEx v8 dataset using `condition-hide=Average` and `greater`.
That differs from the original FLAMES README's 54-specific-tissue command.
The bundled FLAMES `.gsa.out` records the FUMA-style model despite the README
omitting those modifiers. PostGWAS therefore documents both upstream-supported
paths instead of enforcing one as a universal FLAMES requirement. To reproduce
FUMA or the bundled FLAMES result, add `--covariate-model
condition-hide=Average --covariate-direction greater` and use the corresponding
30-general-tissue or 54-specific-tissue covariate table intentionally.

Use `--dry-run` to perform configuration, input, annotation-resource, model,
runtime, and command validation without invoking either upstream analysis stage.

## Parameters

Canonical defaults are defined only in
`config/defaults/modules/flames.yaml` and its typed schema. Omitted CLI options
do not introduce independent defaults. The CLI can override input paths, genome
build, model directory, and explicit `api` or `local` VEP/CADD modes. It does
not expose model distance or weighting controls because changing the published
calibration would require a separately validated model.

## Processing steps

1. Resolve and validate YAML plus explicit CLI overrides once.
2. Verify cross-module genome builds, runtime dependencies, annotation-bundle
   structure, model, and feature manifest.
3. Validate the indexed credible sets and exact MAGMA/PoPS gene compatibility.
4. Write the resolved configuration and a staged index with explicit annotation
   destinations.
5. Run upstream annotation, then verify every expected feature is finite and
   every annotated gene exists in both scientific input sets.
6. Run upstream scoring; verify schemas, ranges, binary indicators, per-locus
   score normalization, and exact agreement between prioritized and causal rows.
7. Rewrite staged annotation paths to published paths, publish validated files,
   and write the completion manifest last.

## Outputs

Default YAML publishes:

- `results/<dataset>_FLAMES_scores.raw`: all scored locus-gene rows.
- `results/<dataset>_FLAMES_scores.pred`: prioritized rows; a header-only file is
  valid when FLAMES prioritizes no gene.
- `annotations/<dataset>/`: one validated annotation table per indexed locus.
- `inputs/<dataset>_flames_index.tsv`: the final index pointing to published
  annotation files.
- `run_metadata/resolved_config.yaml`: effective configuration.
- `run_metadata/<dataset>_flames_completion.yaml`: input/output fingerprints and
  completion status.
- `logs/<dataset>_flames.log`: commands, decisions, counts, failures, and outputs.

Files without a valid completion manifest are incomplete and are not resumed.

## QC and logs

Review cumulative PIP per credible set; chromosome and build declarations;
MAGMA, PoPS, and shared-gene counts; annotation and feature counts; features
that are zero across every locus; scored and prioritized gene counts; resolved
API/local annotation modes; and the completion manifest. A feature that is zero
for every locus is recorded for review because it may indicate absent evidence
or an incompatible annotation resource, but is not automatically scientifically
invalid.

## Interpretation

`FLAMES_scaled` is normalized within a locus and supports ranking genes in that
locus. It is not a posterior probability that a gene is causal. The reported
estimated cumulative precision reflects calibration on the upstream benchmark,
so transferability depends on similarity to that calibration setting and on the
completeness and compatibility of the supplied evidence.

## Common problems

Execution stops before analysis for missing resources, malformed or incomplete
credible sets, build disagreement, duplicate or invalid gene identifiers,
non-finite statistics, no exact MAGMA/PoPS gene overlap, absent dependencies, or
incompatible local VEP/CADD resources. Existing outputs require `--overwrite`
unless `--resume` can validate the complete manifest and all fingerprints.

## Limitations

API mode uses the live Ensembl VEP and CADD services selected by the unchanged
upstream program. Service content and availability may change, so API-mode runs
are not fully reference-pinned even though PostGWAS records the selected mode.
Use build-compatible local resources when reproducible annotation snapshots are
required. PostGWAS validates the wrapper contract and outputs but does not alter
or re-derive the published FLAMES model.

## Scientific references

- [Schipper et al., FLAMES, Nature Genetics (2025)](https://doi.org/10.1038/s41588-025-02084-7)
- [Official FLAMES implementation](https://github.com/Marijn-Schipper/FLAMES)
- [Original FLAMES MAGMA command](https://github.com/Marijn-Schipper/FLAMES/blob/d411689cfe617a5e8bc719e457dfba05e3f7bd3e/README.md#2-run-magma-tissue-type-analysis-using-your-magma-z-scores-on-the-preformatted-gtex-tissue-expression-file)
- [Bundled FLAMES MAGMA result metadata](https://github.com/Marijn-Schipper/FLAMES/blob/d411689cfe617a5e8bc719e457dfba05e3f7bd3e/example_data/magma_exp_gtex_v8_ts_avg_log2TPM.txt.gsa.out#L1-L4)
- [FLAMES annotation and model resources](https://zenodo.org/records/12635505)
- [FUMA FLAMES wrapper and its MAGMA/PoPS inputs](https://github.com/vufuma/FUMA-webapp/blob/0c0259b7ed5e6d15369e978530b7ba56b5bd0437/scripts/flames/run_flames.py#L198-L244)
- [FUMA MAGMA tissue-expression model](https://github.com/vufuma/FUMA-webapp/blob/0c0259b7ed5e6d15369e978530b7ba56b5bd0437/scripts/magma/magma.py#L118-L126)
- [Ensembl VEP REST documentation](https://rest.ensembl.org/documentation/info/vep_region_get)
- [CADD web service and data releases](https://cadd.gs.washington.edu/)
