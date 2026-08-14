# K-POPS

PostGWAS integrates [K-POPS v1.0.0](https://github.com/JasonTan-code/k-pops/tree/8acd49ed8c96565b17c2997420f514f24b65e096), a kernel-based gene-prioritisation method that fits kernel ridge regression to MAGMA gene statistics and reports supporting contributor genes. The Docker build pins upstream commit `8acd49ed8c96565b17c2997420f514f24b65e096` and runs it with the main `postgwas` conda environment.

## Method contract

K-POPS does not infer genome build or gene-identifier compatibility. Configure one build shared by:

- the MAGMA gene association;
- the K-POPS gene annotation;
- the genes represented by the binary kernel matrix.

The annotation must contain the configured Ensembl ID, gene name, chromosome, and TSS columns. Ensembl IDs must be unique; duplicate gene symbols are permitted because they occur in the official upstream annotation. A gene-symbol anchor is rejected when it maps to multiple Ensembl IDs, avoiding the upstream dictionary conversion's otherwise ambiguous choice. PostGWAS verifies kernel dimensions, finite values, and symmetry in bounded row chunks before analysis. It also rejects anchor genes that upstream K-POPS would otherwise discard and verifies that every requested fit has enough training genes for the configured contributor count.

`gene_universe_policy` controls MAGMA genes that are absent from the K-POPS annotation or kernel. The default `intersect` policy retains only genes present in all three inputs, preserves their original MAGMA order, subsets `.genes.out` and `.genes.raw` together in temporary run storage, warns with retained and excluded counts, and publishes a gene-level TSV plus YAML compatibility report. Original MAGMA files are never changed. The temporary `.genes.raw` subset is used only for the technical MAGMA covariates read by K-POPS and is deleted after the run. Set the policy to `strict` to stop before fitting whenever any MAGMA gene is incompatible. A chromosome disagreement for a shared gene always fails because an identifier intersection cannot repair inconsistent genomic annotation.

This default makes upstream K-POPS's implicit target selection explicit and auditable. The pinned upstream implementation intersects MAGMA targets with annotation chromosomes before it indexes the kernel; PostGWAS additionally validates the retained universe and refuses to run when it falls below `minimum_gene_count`.

Upstream K-POPS writes its gene identifier as an unnamed table index. PostGWAS performs one explicit publication transformation: that index is named with `input_schema.published_gene_id_column` (default `ENSGID`). Scores and gene order are unchanged. The transformation is recorded by the resolved configuration and output manifest.

The upstream prediction table is parsed with its literal tab delimiter so consecutive tabs remain explicit missing fields. This is distinct from the whitespace-delimited annotation and MAGMA inputs; treating predictions as generic whitespace would collapse missing support fields and corrupt row boundaries.

After validating the published predictions, PostGWAS prints and logs an analysis summary containing the finite-score coverage, MAGMA target count, configured training design, top-ranked genes, interpretation guidance, and chromosome-level coverage warnings. `reporting.top_gene_count` and `reporting.score_decimal_places` control presentation only and do not change the K-POPS fit or invalidate resumable analysis outputs. Genes without finite scores remain in the published table but are excluded from the displayed ranking. K-POPS coefficients describe kernel-regression training genes and are not reported as selected biological features.

The upstream program defines its HLA interval internally. PostGWAS exposes only the upstream keep/remove switches; it does not redefine or reinterpret that interval.

Direct execution reports four validated progress stages: input and reference validation, compatible MAGMA gene-universe resolution, K-POPS model fitting and scoring, and final output validation and publication. The model-fitting stage remains active with elapsed time while upstream K-POPS runs; PostGWAS does not invent chromosome-level completion percentages that the upstream program does not expose.

## Inputs

- `magma_association_prefix`: prefix for `.genes.out` and `.genes.raw`.
- `gene_annotation_file`: K-POPS annotation containing the configured columns.
- `kernel_matrix_prefix`: prefix for the float32 `.bin` matrix and `.genes` order file.
- `gene_universe_policy`: `intersect` (default) or `strict` handling for incompatible MAGMA target genes.
- `script_path`: installed `k-pops.py` command name, normally resolved automatically from the active environment; `--kpops-script` may override it for a nonstandard installation.

Use upstream `prepare_kernel.py` to construct a kernel from PoPS feature matrices. Kernel type, standardisation, RBF gamma, and polynomial degree are resource-preparation decisions and must be recorded when the kernel is created; PostGWAS does not silently rebuild an analysis kernel.

## Key command controls

- `--training-chromosomes loco` fits leave-one-chromosome-out models; `all` uses all chromosomes for training and prediction; explicit annotation chromosome labels train one model on that subset and predict the remaining chromosomes.
- `--use-magma-covariates` projects gene-size, gene-density, and inverse-MAC covariates derived from MAGMA `.genes.raw` out of the target Z statistics before fitting. `--no-use-magma-covariates` disables that adjustment.
- `--top-contributor-gene-count` controls how many highest positive contribution genes appear in each prediction's explanation columns; it does not alter the score.
- `--anchor-genes` accepts space-separated Ensembl IDs or gene symbols according to `--anchor-gene-type` and enables `Anchor_Score` evidence. Ambiguous gene symbols are rejected.
- `--save-attribution-files` writes the dense contribution matrix and row/column gene-order files. Its storage requirement grows quadratically with the number of genes.
- `--kpops-device` selects the PyTorch device: `cpu`, `cuda` for a supported NVIDIA GPU, or `mps` for supported Apple Silicon.

## Direct mode

```console
postgwas kpops \
  --magma-association-prefix magma/02_intermediates/positional/native_outputs/STUDY_magma_35up_10down \
  --kpops-gene-annotation-file /path/to/postgwas-resources/pops/GRCh37_gene_annot_jun10.txt \
  --kernel-matrix-prefix /path/to/postgwas-resources/kpops/kernels/GRCh37/pops_features_standardized_linear \
  --kpops-genome-build GRCh37 \
  --dataset-id STUDY \
  --output-directory results
```

The following command translates every model option in the upstream v1.0.0 ApoB README example into the PostGWAS CLI. `--seed 42` preserves the upstream program's default seed, while `--gene-universe-policy strict` prevents PostGWAS's auditable intersection policy from changing the supplied ApoB gene universe. The upstream repository does not publish the `ApoB.genes.out`, `ApoB.genes.raw`, or prepared `test/kernel_linear` inputs, so `/path/to/ApoB` and `/path/to/ApoB/kernel_linear` must point to copies prepared from the original study inputs. CUDA also requires a supported NVIDIA GPU.

```console
postgwas kpops \
  --magma-association-prefix /path/to/ApoB \
  --kpops-gene-annotation-file /path/to/postgwas-resources/pops/GRCh37_gene_annot_jun10.txt \
  --kernel-matrix-prefix /path/to/ApoB/kernel_linear \
  --kpops-genome-build GRCh37 \
  --gene-universe-policy strict \
  --use-magma-covariates \
  --training-chromosomes loco \
  --kpops-device cuda \
  --top-contributor-gene-count 5 \
  --anchor-genes PCSK9 ANGPTL3 APOB ABCG5 ALB NPC1L1 GIGYF1 JAK2 A1CF PDE3B APOC3 CETP ASGR1 LDLR ZNF234 ZNF229 NECTIN2 RRBP1 \
  --anchor-gene-type NAME \
  --seed 42 \
  --kpops-verbose \
  --dataset-id ApoB \
  --output-directory results/ApoB
```

## Pipeline mode

```console
postgwas config export --pipeline kpops --style full --output kpops_pipeline.yaml

postgwas pipeline \
  --modules kpops \
  --vcf study_GRCh37.vcf.gz \
  --magma-ld-reference /path/to/postgwas-resources/magma/functional_mapping/base/ld_reference/g1000_eur/g1000_eur \
  --gene-location-file /path/to/postgwas-resources/kpops/software/8acd49ed8c96565b17c2997420f514f24b65e096/data/Ensembl.hg19.gene.loc \
  --kpops-gene-annotation-file /path/to/postgwas-resources/pops/GRCh37_gene_annot_jun10.txt \
  --kernel-matrix-prefix /path/to/postgwas-resources/kpops/kernels/GRCh37/pops_features_standardized_linear \
  --kpops-genome-build GRCh37 \
  --dataset-id STUDY \
  --output-directory results \
  --run-config kpops_pipeline.yaml
```

Pipeline mode runs formatter and MAGMA first, then supplies the validated MAGMA prefix to K-POPS. The displayed `Ensembl.hg19.gene.loc` is the file pinned with K-POPS v1.0.0; it uses Ensembl identifiers and the same GRCh37 gene universe as the installed K-POPS annotation. An NCBI/Entrez MAGMA location file is not compatible with this kernel. Installation provides the pinned `k-pops.py` command. The annotation, kernel prefix, and declared build must be supplied explicitly through the corresponding CLI options or a run configuration; PostGWAS has no machine-specific data-resource defaults.

## Outputs

K-POPS results are staged and published only after validation:

- `{dataset_id}_kpops.preds`: canonical `ENSGID` and `PoPS_Score` predictions plus explanation columns;
- `{dataset_id}_kpops.coefs`: fitted coefficients;
- optional attribution matrix and row/column gene files;
- with `gene_universe_policy: intersect`, a gene-level compatibility TSV and YAML report;
- resolved configuration, completion manifest, and canonical log.

## Sources

- [K-POPS v1.0.0 README](https://github.com/JasonTan-code/k-pops/tree/8acd49ed8c96565b17c2997420f514f24b65e096)
- [Pinned K-POPS implementation](https://github.com/JasonTan-code/k-pops/blob/8acd49ed8c96565b17c2997420f514f24b65e096/k-pops.py)
