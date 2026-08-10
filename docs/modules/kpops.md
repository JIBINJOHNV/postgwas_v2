# K-POPS

PostGWAS integrates [K-POPS v1.0.0](https://github.com/JasonTan-code/k-pops/tree/8acd49ed8c96565b17c2997420f514f24b65e096), a kernel-based gene-prioritisation method that fits kernel ridge regression to MAGMA gene statistics and reports supporting contributor genes. The Docker build pins upstream commit `8acd49ed8c96565b17c2997420f514f24b65e096` and runs it with the main `postgwas` conda environment.

## Scientific contract

K-POPS does not infer genome build or gene-identifier compatibility. Configure one build shared by:

- the MAGMA gene association;
- the K-POPS gene annotation;
- the genes represented by the binary kernel matrix.

The annotation must contain the configured Ensembl ID, gene name, chromosome, and TSS columns. Ensembl IDs must be unique; duplicate gene symbols are permitted because they occur in the official upstream annotation. A gene-symbol anchor is rejected when it maps to multiple Ensembl IDs, avoiding the upstream dictionary conversion's otherwise ambiguous choice. MAGMA gene IDs must occur exactly in both the annotation and kernel. PostGWAS verifies kernel dimensions, finite values, and symmetry in bounded row chunks before analysis. It also rejects anchor genes that upstream K-POPS would otherwise discard and verifies that every requested fit has enough training genes for the configured contributor count.

Upstream K-POPS writes its gene identifier as an unnamed table index. PostGWAS performs one explicit publication transformation: that index is named with `input_schema.published_gene_id_column` (default `ENSGID`). Scores and gene order are unchanged. The transformation is recorded by the resolved configuration and output manifest.

The upstream prediction table is parsed with its literal tab delimiter so consecutive tabs remain explicit missing fields. This is distinct from the whitespace-delimited annotation and MAGMA inputs; treating predictions as generic whitespace would collapse missing support fields and corrupt row boundaries.

After validating the published predictions, PostGWAS prints and logs a scientific summary containing the finite-score coverage, MAGMA target count, configured training design, top-ranked genes, interpretation guidance, and chromosome-level coverage warnings. `reporting.top_gene_count` and `reporting.score_decimal_places` control presentation only and do not change the K-POPS fit or invalidate resumable analysis outputs. Genes without finite scores remain in the published table but are excluded from the displayed ranking. K-POPS coefficients describe kernel-regression training genes and are not reported as selected biological features.

The upstream program defines its HLA interval internally. PostGWAS exposes only the upstream keep/remove switches; it does not redefine or reinterpret that interval.

## Inputs

- `magma_association_prefix`: prefix for `.genes.out` and `.genes.raw`.
- `gene_annotation_file`: K-POPS annotation containing the configured columns.
- `kernel_matrix_prefix`: prefix for the float32 `.bin` matrix and `.genes` order file.
- `script_path`: installed `k-pops.py` command name, normally resolved automatically from the active environment; `--kpops-script` may override it for a nonstandard installation.

Use upstream `prepare_kernel.py` to construct a kernel from PoPS feature matrices. Kernel type, standardisation, RBF gamma, and polynomial degree are resource-preparation decisions and must be recorded when the kernel is created; PostGWAS does not silently rebuild an analysis kernel.

## Direct mode

```console
postgwas kpops \
  --magma-association-prefix magma/results/STUDY_magma_35up_10down \
  --kpops-gene-annotation-file reference/kpops/GRCh37_gene_annot.tsv \
  --kernel-matrix-prefix reference/kpops/kernel_linear \
  --kpops-genome-build GRCh37 \
  --dataset-id STUDY \
  --output-directory results
```

## Pipeline mode

```console
postgwas config export --pipeline kpops --style full --output kpops_pipeline.yaml

postgwas pipeline \
  --modules kpops \
  --vcf study_GRCh37.vcf.gz \
  --magma-ld-reference reference/1000G_EUR \
  --gene-location-file reference/NCBI37.3.gene.loc \
  --kpops-gene-annotation-file reference/kpops/GRCh37_gene_annot.tsv \
  --kernel-matrix-prefix reference/kpops/kernel_linear \
  --kpops-genome-build GRCh37 \
  --dataset-id STUDY \
  --output-directory results \
  --run-config kpops_pipeline.yaml
```

Pipeline mode runs formatter and MAGMA first, then supplies the validated MAGMA prefix to K-POPS. Installation provides the pinned `k-pops.py` command. The annotation, kernel prefix, and declared build must be supplied explicitly through the corresponding CLI options or a run configuration; PostGWAS has no machine-specific data-resource defaults.

## Outputs

K-POPS results are staged and published only after validation:

- `{dataset_id}_kpops.preds`: canonical `ENSGID` and `PoPS_Score` predictions plus explanation columns;
- `{dataset_id}_kpops.coefs`: fitted coefficients;
- optional attribution matrix and row/column gene files;
- resolved configuration, completion manifest, and canonical log.

## Sources

- [K-POPS v1.0.0 README](https://github.com/JasonTan-code/k-pops/tree/8acd49ed8c96565b17c2997420f514f24b65e096)
- [Pinned K-POPS implementation](https://github.com/JasonTan-code/k-pops/blob/8acd49ed8c96565b17c2997420f514f24b65e096/k-pops.py)
