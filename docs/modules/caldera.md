# CALDERA

PostGWAS integrates [CALDERA](https://github.com/kheilbron/caldera/tree/81a8a0308741ae986660f711bbf6a7abbd3bca19), which prioritises causal genes within GWAS loci using PoPS score, distance, coding posterior inclusion probability, and a trained logistic model. The Docker build pins upstream commit `81a8a0308741ae986660f711bbf6a7abbd3bca19` in the main `postgwas` conda environment.

## Scientific contract

CALDERA consumes:

- PoPS predictions with `ENSGID` and `PoPS_Score`;
- one credible-set variant per row with `locus`, `chr`, `bp`, and `pip`;
- an explicit GRCh37 or GRCh38 declaration.

PostGWAS rejects missing required values, non-integer chromosome or position labels, PIPs outside `[0, 1]`, loci spanning multiple chromosomes, duplicate PoPS genes, non-finite PoPS scores, and credible sets below 95% cumulative PIP. This matches the pinned upstream implementation: totals exactly at 0.95 remain valid, while totals above the boundary are retained through the first row that makes cumulative PIP greater than 0.95. PostGWAS records how many supplied rows upstream will retain and trim.

The `minimum_credible_set_coverage` configuration key is schema-locked to `0.95`. It records the upstream protocol invariant rather than offering an ineffective override: the pinned CALDERA R function fixes this boundary internally and exposes no cutoff argument.

The upstream coding-variant position resource is GRCh37. Therefore:

- GRCh37 credible sets may be annotated by chromosome and position;
- GRCh38 direct-mode inputs must supply the configured `c_gene` or `rsid` column;
- current pipeline conversion is supported for GRCh37 only because the shared fine-mapping interchange does not retain `c_gene` or `rsid`.

PostGWAS never lifts coordinates or treats K-POPS scores as PoPS scores. CALDERA pipeline mode depends specifically on the PoPS module because the trained CALDERA model was defined for PoPS evidence.

Installation places the pinned CALDERA repository under the active `postgwas` environment, and PostGWAS resolves it relative to that environment's Python. The R adapter is resolved from the installed PostGWAS package. `--caldera-repository` and `--caldera-adapter-script` remain optional overrides for nonstandard development installations.

## Direct mode

```console
postgwas caldera \
  --pops-file pops/STUDY_pops.preds \
  --credible-set-file finemap/STUDY_credible_sets_GRCh37.tsv \
  --caldera-genome-build GRCh37 \
  --dataset-id STUDY \
  --output-directory results
```

For GRCh38, include `c_gene` or `rsid` in the credible-set table. Column names are controlled by `modules.caldera.input_schema`.

## Pipeline mode

```console
postgwas config export --pipeline caldera --style full --output caldera_pipeline.yaml

postgwas pipeline \
  --modules caldera \
  --vcf study_GRCh37.vcf.gz \
  --genome-build GRCh37 \
  --ld-region-dir reference/ld_blocks \
  --ld-block-populations EUR \
  --ld-folder reference/ld \
  --population EUR \
  --magma-ld-reference reference/1000G_EUR \
  --gene-location-file reference/NCBI37.3.gene.loc \
  --feature-matrix-prefix reference/pops/features_munged/pops_features \
  --feature-matrix-chunks 116 \
  --pops-gene-location-file reference/pops/GRCh37_gene_annot.tsv \
  --pops-genome-build GRCh37 \
  --finemap-method susie \
  --finemap-ld-reference reference/1000G_EUR \
  --caldera-genome-build GRCh37 \
  --dataset-id STUDY \
  --output-directory results \
  --run-config caldera_pipeline.yaml
```

The plan runs MAGMA and PoPS plus fine-mapping prerequisites. SuSiE and FINEMAP already publish validated 95% credible sets in the common FLAMES interchange. CALDERA converts that artifact to `locus/chr/bp/pip`, using the configured index filenames, column names, variant-ID pattern, and locus separator. The converted table is retained with the result for provenance.

## Outputs

- `results/{dataset_id}_caldera.tsv`: gene probabilities and component evidence;
- `inputs/{dataset_id}_caldera_credible_sets.tsv`: pipeline conversion, when applicable;
- resolved configuration, completion manifest, and canonical log.

CALDERA `caldera` probabilities are validated to sum to one per locus. The `multi` column is the upstream unnormalised probability.

## Sources

- [CALDERA README and input/output specification](https://github.com/kheilbron/caldera/tree/81a8a0308741ae986660f711bbf6a7abbd3bca19)
- [Pinned CALDERA implementation](https://github.com/kheilbron/caldera/blob/81a8a0308741ae986660f711bbf6a7abbd3bca19/z_caldera.R)
- [CALDERA preprint](https://www.medrxiv.org/content/10.1101/2024.07.26.24311057v2)
