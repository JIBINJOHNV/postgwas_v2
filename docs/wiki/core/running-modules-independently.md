# Running Modules Independently

Direct mode runs a named analysis using that module's inputs. It does not ask the
global pipeline planner to schedule all upstream modules. You supply the
required VCF, formatted tables, or upstream results, plus the relevant references.
The module still performs its own input and compatibility checks; direct mode is
not an unvalidated shortcut.

Use [Pipeline Workflow](pipeline-workflow.md) when PostGWAS should arrange and
pass the registered upstream artifacts. Neither mode can make incompatible
studies and references suitable for a method.

## Standalone preparation

`harmonisation` is standalone. Use it to convert mapped raw summary-statistics
columns into the PostGWAS-harmonised GWAS-VCF and QC evidence. If you already
have a compatible PostGWAS-produced VCF, you do not repeat preparation merely
to switch between direct and pipeline execution. Arbitrary third-party VCFs
do not meet the [study-VCF provenance boundary](input-output-contracts.md#postgwas-origin-check-for-study-vcfs).

## Downstream standalone commands

Every analysis module has a standalone command: `sumstat_filter`, `formatter`,
`imputation`, `annot_ldblock`, `ld_clump`, `manhattan`, `qc`, `heritability`,
`finemap`, `magma`, `gcta_cojo`, `gcta_gene`, `magmacovar`, `single_cell`,
`pops`, `kpops`, `caldera`, `flames`, `mixer`, and `pathway_enrichment`. Their
exact installed names are shown by `postgwas --help`.

Before running one directly:

1. Read its `--help` output from the installed checkout.
2. Confirm the required upstream artifact exists and was produced with matching
   build, ancestry, alleles, and sample-size definitions.
3. Review the resolved settings. Use an exported run configuration where the
   command exposes `--run-config`; do not assume every direct parser accepts it.
4. Verify external executables and reference files.
5. Use a new output directory or dataset identifier to avoid confusing outputs
   from different runs.
6. Retain the command, configuration, logs, and QC summaries.

## Example: MAGMA with and without pipeline orchestration

These two routes use the GRCh37/EUR/NCBI37.3 example from
[Quick Start](../getting-started/quick-start.md), with the same MAGMA settings
and reference files. For the explicit formatter example below, the reference
BIM uses rsIDs. A reference using supported coordinate/allele identifiers needs
the corresponding formatter selection instead; review the identifier audit.

### Direct: prepare the tables, then run MAGMA

```console
postgwas formatter \
  --vcf results/harmonisation/STUDY/harmonisation/STUDY_GRCh37_merged.vcf.gz \
  --format magma \
  --variant-id-type rsid \
  --dataset-id STUDY \
  --output-directory results/formatted

postgwas magma \
  --snp-location-file results/formatted/STUDY_magma_snp_loc.tsv \
  --p-value-file results/formatted/STUDY_magma_p_values.tsv \
  --magma-ld-reference reference/GRCh37_EUR_reference \
  --gene-location-file reference/NCBI37.3.gene.loc \
  --resolve-variants-to-reference \
  --dataset-id STUDY \
  --output-directory results/magma_direct
```

### Pipeline: request MAGMA from the VCF

```console
postgwas pipeline \
  --modules magma \
  --vcf results/harmonisation/STUDY/harmonisation/STUDY_GRCh37_merged.vcf.gz \
  --magma-ld-reference reference/GRCh37_EUR_reference \
  --gene-location-file reference/NCBI37.3.gene.loc \
  --resolve-variants-to-reference \
  --dataset-id STUDY \
  --output-directory results/magma_pipeline
```

The pipeline validates the BIM identifier convention and configures the MAGMA
formatter target accordingly. It then runs `formatter → magma` and passes the
created tables automatically. Both commands still need the external LD and gene
references; both apply the explicitly requested exact-ID reference filtering.
Compare resolved settings and retention counts, not directory names, when
comparing runs. See [MAGMA](../modules/magma.md) for analysis policies, outputs,
and optional competitive gene-set analysis.

## Method-dependent input relationships

- MAGMA, GCTA gene tests, GCTA-COJO, LDSC, MiXeR, PRED-LD imputation, and
  fine-mapping consume their own formatter contracts, not an arbitrary table.
- Standard LD clumping uses the harmonised VCF and a prepared pairwise-LD
  reference; it does not require LD-block annotation. Region pruning does
  require that annotation. Direct `ld_clump` with `cojo-slct` can prepare its
  GCTA table internally: direct mode does not prohibit a module's own input
  preparation.
- Direct fine-mapping requires a locus file plus the selected engine's
  formatted table and LD reference. The pipeline instead obtains its locus
  artifact from standard clumping.
- MAGMAcovar, PoPS, and K-POPS consume MAGMA results. The single-cell
  MAGMA/scDRS routes also use MAGMA in the pipeline, whereas `ldsc_celltype`
  depends on the LDSC formatter route instead.
- CALDERA and FLAMES consume several compatible upstream results; follow their
  method-specific guides rather than treating them as independent VCF readers.

MAGMA cell typing consumes an exact MAGMA `.genes.raw` file and an aggregated
single-cell expression covariate matrix in direct mode. Pipeline mode supplies
the gene results through the MAGMA dependency, but the biological covariate
matrix remains an explicit user input.

Pathway enrichment and harmonisation are standalone-only.

## Configuration and execution exceptions

Most direct analysis commands accept `--run-config`, but direct `manhattan`
and `pathway_enrichment` currently do not. Pipeline Manhattan does accept the
pipeline run configuration. Pathway enrichment takes its gene list and provider
credentials directly, runs a fixed provider workflow, and cannot be selected as
a pipeline target. Its documented provider configuration is not a provider
selection interface.

Use the exact public command's help and [Configuration](configuration.md) to
identify which settings are exposed. Preserve the full command and any supplied
configuration with the [output record](../reference/output-structure.md).
