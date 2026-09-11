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

Direct and pipeline execution report six validated K-POPS stages. PostGWAS first validates the two MAGMA gene-association input files, then validates the gene annotation, then the kernel matrix, gene-order file, and installed K-POPS script. Only after those independent checks pass does it validate the cross-resource gene universe, run model fitting, and validate and publish the outputs. Each completed validation stage prints the relevant paths, row or gene counts, declared build, and structural result; the same evidence is retained in the canonical log.

During stage 5, a separate `K-POPS model-fitting progress` counter reports observed upstream model completions when detailed progress is enabled. PostGWAS enables the pinned program's verbose messages internally and captures them rather than printing duplicate native output. Each upstream `Computing PoPS scores.` event occurs after the corresponding coefficient solve and advances the measured counter. LOCO analysis uses the exact number of distinct chromosomes in the validated K-POPS annotation as its total; `all` and explicitly selected training chromosomes each fit one model. This measured counter is distinct from the six-stage PostGWAS bar: for example, `5/6` identifies the active PostGWAS stage, while `4/6 67%` means four stages have completed validation. The measured counter reaches 100% only after the K-POPS process exits successfully and its declared files pass the checked-command output contract. Full scientific result validation remains stage 6. A failed run retains its last observed count below 100%. Enabling verbose capture changes logging only, not the K-POPS scientific arguments or calculations.

## Input requirements

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

The following command needs no run YAML. Replace paths and uppercase
source-metadata placeholders with the exact reference manifest values; the
source URL must identify the actual gene-location resource over HTTPS. The
Ensembl declaration describes column one; it does not convert identifiers.

```console
postgwas pipeline \
  --modules kpops \
  --vcf study_GRCh37.vcf.gz \
  --magma-ld-reference /path/to/postgwas-resources/magma/functional_mapping/base/ld_reference/g1000_eur/g1000_eur \
  --gene-location-file reference/PoPS_GRCh37_strand_aware.loc \
  --magma-positional-gene-id-type ensembl \
  --magma-positional-source-name GENE_LOCATION_SOURCE \
  --magma-positional-source-version SOURCE_RELEASE \
  --magma-positional-source-url https://example.org/GENE_LOCATION_SOURCE_RECORD \
  --magma-positional-context GENE_LOCATION_CONTEXT \
  --kpops-gene-annotation-file /path/to/postgwas-resources/pops/GRCh37_gene_annot_jun10.txt \
  --kernel-matrix-prefix /path/to/postgwas-resources/kpops/kernels/GRCh37/pops_features_standardized_linear \
  --kpops-genome-build GRCh37 \
  --dataset-id STUDY \
  --output-directory results
```

To keep the same settings in YAML instead, export
`postgwas config export --pipeline kpops --style full --output kpops_pipeline.yaml`,
edit the required inputs and positional mapping declarations, and provide it
with `--run-config`. An unedited export does not supply study-specific resources.

Pipeline mode runs formatter and MAGMA first, then supplies the validated MAGMA prefix to K-POPS. The prepared gene-location reference must use the same GRCh37 Ensembl gene universe as the K-POPS annotation, with columns `gene_id chromosome start end strand` and an optional symbol. PostGWAS requires strand to be `+` or `-`; the upstream K-POPS `Ensembl.hg19.gene.loc` instead contains a numeric TSS in column five and is not directly compatible. An explicitly audited preparation may derive `+` when the source TSS equals START and `-` when it equals END, but must reject ambiguous or non-endpoint TSS values and preserve all source genes and coordinates. An NCBI/Entrez location file is not compatible with this kernel. Installation provides the pinned `k-pops.py` command. The annotation, kernel prefix, and declared build must be supplied explicitly through the corresponding CLI options or a run configuration; PostGWAS has no machine-specific data-resource defaults.

Pipeline-only overrides use `--kpops-gene-universe-policy`, `--kpops-training-chromosomes`, and `--kpops-use-magma-covariates` / `--no-kpops-use-magma-covariates`. The corresponding direct-mode names remain unchanged. This keeps K-POPS settings independent when PoPS is selected in the same pipeline: unqualified PoPS flags do not override K-POPS, and omitted options retain each module's canonical YAML values. No configuration keys or model defaults change.

## Outputs

K-POPS results are staged and published only after validation:

- `{dataset_id}_kpops.preds`: canonical `ENSGID` and `PoPS_Score` predictions plus explanation columns;
- `{dataset_id}_kpops.coefs`: fitted coefficients;
- optional attribution matrix and row/column gene files;
- with `gene_universe_policy: intersect`, a gene-level compatibility TSV and YAML report;
- resolved configuration, completion manifest, and canonical log.

## Interpretation

Use `PoPS_Score` to rank genes under the chosen kernel and training design; it
is neither a p-value nor a calibrated causal probability. LOCO predictions
exclude the prediction chromosome from model fitting; `all` and selected
training-chromosome designs answer different prediction questions and should
not be presented as interchangeable validation. Review finite-score coverage
and the gene-universe compatibility report before comparing rankings.
Contributor genes explain the fitted kernel prediction, not a demonstrated
regulatory interaction or a set of selected biological features.

## Common problems

Check exact gene identifiers, shared-gene chromosome agreement, kernel row order
and dimensions, declared build, and adequate training genes. A symbol-based
MAGMA file or a numerically encoded TSS in the gene-location strand column is
not repaired by adding an Ensembl metadata label. Follow the explicit resource
contract above; PostGWAS does not offer a general gene-location conversion CLI.

## Sources

- [K-POPS v1.0.0 README](https://github.com/JasonTan-code/k-pops/tree/8acd49ed8c96565b17c2997420f514f24b65e096)
- [Pinned K-POPS implementation](https://github.com/JasonTan-code/k-pops/blob/8acd49ed8c96565b17c2997420f514f24b65e096/k-pops.py)
- [MAGMA documentation and strand-aware annotation format](https://cncr.nl/research/magma/)
