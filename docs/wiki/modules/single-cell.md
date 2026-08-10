# Single-Cell Integration

## Purpose

The `single_cell` module integrates gene-level GWAS evidence with single-cell
expression. It currently provides the FUMA-compatible MAGMA cell-type base
analysis and scDRS cell scoring. Select either method, or both in one direct
run, with `--tools magma_celltype scdrs`.

## What the analysis does

The pipeline reuses PostGWAS MAGMA to convert variant-level GWAS evidence into
calibrated gene-association results. It then reuses the MAGMAcovar service to
regress the MAGMA gene statistic on each cell-type average-expression property
while conditioning on the configured across-cell-type average-expression
property. The alternative hypothesis is prespecified as a positive association.

For cell type \(c\), the FUMA base model is conceptually
\(Z = \beta_0 + E_c\beta_E + A\beta_A + B\beta_B + \epsilon\), where \(Z\)
is the MAGMA gene statistic, \(E_c\) is average expression in cell type \(c\),
\(A\) is average expression across cell types, and \(B\) represents MAGMA's
technical gene covariates. The implemented engine call is equivalent to
`--model condition-hide=Average direction-covar=greater` with the property name
resolved from configuration. FUMA writes the general modifier as
`direction=greater`; the reused MAGMAcovar service uses MAGMA's documented
covariate-specific spelling, `direction-covar=greater`, for the same positive
one-sided covariate hypothesis.

scDRS answers a different question. It selects a disease gene set from GWAS
gene statistics, then scores every cell by comparing expression of those genes
with matched control gene sets. Native output includes raw and normalized cell
scores, Monte Carlo and pooled p-values, and optional cell-group association,
within-group heterogeneity, continuous-annotation correlation, and gene-score
correlation results. PostGWAS preserves those native values; it does not call
the per-cell score a cell-type association or silently adjust its raw p-values.

## When to use it

Use `magma_celltype` when the question is whether cell-type average-expression
properties predict MAGMA gene association. Use `scdrs` when the question is
which individual cells, annotated groups, or cell states preferentially express
a GWAS-derived disease gene set. Use pipeline mode when starting from a
harmonised GWAS-VCF. Use direct mode with exact method-specific inputs.

## Input requirements

For `magma_celltype`:

- Pipeline mode requires a harmonised GWAS-VCF and the references required by
  the existing MAGMA module; direct mode requires compatible `.genes.raw`.
- Both modes require a whitespace-delimited expression covariate matrix with
  gene IDs, one column per cell type, and the configured average-expression
  column (packaged name `Average`).
- Gene identifiers must match the MAGMA result.

For `scdrs`:

- Direct file mode requires an H5AD file plus a native two-column scDRS `.gs`
  file. The header is `TRAIT` and `GENESET`; genes may have `gene:weight` values.
- Pipeline mode requires an H5AD file and, by default, a pinned one-to-one
  crosswalk from the primary MAGMA gene identifiers to `adata.var_names`.
  PostGWAS obtains the headered `.genes.out` from its MAGMA stage and constructs
  the `.gs`; users do not provide a disease-relevance score file.
- The configured expression matrix is exactly `adata.X`. Declare it as
  `raw_counts` or `normalized_log1p`. Raw-count mode rejects fractional values;
  both modes reject non-finite and negative values.
- H5AD cell and gene identifiers must be unique. Requested group and continuous
  annotations must exist in `adata.obs`; continuous annotations must be numeric.
- A crosswalk must use the configured source and target columns and be one-to-one.
  Ambiguous source or target mappings fail. Unmapped source genes follow the
  explicit `unmapped_policy` and are counted in the mapping report.
- After filtering and H5AD overlap, every gene set must contain at least 51
  genes and cover less than 20% of the expression-gene universe under packaged
  settings.
- An optional scDRS `.cov` file must contain numeric covariates keyed by cell ID;
  packaged settings require its cell IDs to match `adata.obs_names` exactly and
  require a `const` column containing ones. If exact matching is disabled, cell
  overlap must still be greater than the configured 75% threshold used by
  upstream scDRS.
- The current validator requires the H5AD and gene-set species to be the same.
  Native cross-species homolog conversion is not exposed until PostGWAS can
  validate the post-conversion effective gene set itself.

The expression covariate matrix and H5AD atlas are processed biological inputs.
The GWAS formatter does not create them, and this module does not perform FASTQ
processing, alignment, study-level cell QC, clustering, or annotation.

## Command

```console
postgwas single_cell --tools magma_celltype --magma-gene-results-file PATH \
  --single-cell-covariates PATH [options]

postgwas single_cell --tools scdrs --scdrs-h5ad-file ATLAS.h5ad \
  --scdrs-gene-set-file TRAIT.gs [options]

postgwas single_cell --tools scdrs --scdrs-gene-set-source magma \
  --scdrs-h5ad-file ATLAS.h5ad \
  --scdrs-magma-gene-results-file STUDY.genes.out \
  --scdrs-gene-id-map entrez_to_symbol.tsv [options]
```

## Minimal example

```console
postgwas single_cell \
  --tools magma_celltype \
  --magma-gene-results-file STUDY.genes.raw \
  --single-cell-covariates atlas_celltype_average.tsv \
  --magma /path/to/magma \
  --dataset-id STUDY \
  --output-directory results
```

## Full example

```console
postgwas pipeline \
  --modules single_cell \
  --tools magma_celltype \
  --vcf study.gwas.vcf.gz \
  --single-cell-covariates atlas_celltype_average.tsv \
  --magma-ld-reference reference/g1000_eur \
  --gene-location-file reference/NCBI37.3.gene.loc \
  --cell-type-correction bonferroni fdr_bh \
  --primary-cell-type-correction bonferroni \
  --cell-type-significance-threshold 0.05 \
  --dataset-id STUDY \
  --output-directory results
```

The pipeline stages are `formatter`, `magma`, and `single_cell`. The formatter
creates MAGMA SNP and p-value tables from the GWAS-VCF; MAGMA produces one
validated primary `.genes.raw`; and the single-cell stage runs and normalizes
the cell-type gene-property analysis.

For pipeline scDRS, the same stages are used, but the final stage consumes
MAGMA's headered gene output and creates a weighted disease gene set first:

```console
postgwas pipeline \
  --modules single_cell \
  --tools scdrs \
  --vcf study.gwas.vcf.gz \
  --scdrs-h5ad-file brain_atlas.h5ad \
  --scdrs-gene-id-map entrez_to_symbol.tsv \
  --scdrs-group-analysis cell_type \
  --magma-ld-reference reference/g1000_eur \
  --gene-location-file reference/NCBI37.3.gene.loc \
  --dataset-id STUDY \
  --output-directory results
```

The packaged crosswalk schema expects `ENTREZID` and `SYMBOL`. Override the
columns and the declared source/target identifier types in YAML when using a
different pinned namespace. Set `scdrs.magma_gene_set.mapping_mode: exact` only
when MAGMA and `adata.var_names` genuinely use the same identifiers.

## Parameters

`modules.single_cell.tools` is the canonical ordered tool selection and is
overridden by `--tools`. The accepted values are `magma_celltype` and `scdrs`.
Each tool validates only its own inputs. MAGMA cell typing retains its normalized
result schema and multiple-testing family. scDRS has separate H5AD, gene-set,
mapping, filtering, control-set, downstream-analysis, version, and native-output
settings under `modules.single_cell.scdrs`.

`magma_celltype.magmacovar_use_case` selects a named definition from
`modules.magmacovar.model_use_cases`. PostGWAS verifies that the selected use
case resolves to the FUMA base model and positive one-sided direction. The
remaining MAGMAcovar missingness, overlap, executable, staging, and validation
settings continue to be owned and logged by the reused MAGMAcovar service.

## Processing steps

1. Resolve and schema-validate single-cell, MAGMA, and MAGMAcovar settings.
2. Validate the covariate matrix, exact `Average` property, cell-type columns,
   numeric values, missingness, property variation, and minimum gene count.
3. In pipeline mode, reuse the calibrated primary `.genes.raw` published by the
   preceding MAGMA stage; chromMAGMA ranking output is rejected.
4. Invoke the existing MAGMAcovar service with the configured tissue-
   specificity use case in an isolated engine directory.
5. Require the MAGMA COVAR results to match the complete set of input cell-type
   properties exactly. The conditioning average is not part of the tested
   cell-type family.
6. Apply every configured correction across all tested cell types in that
   dataset without filtering or reordering the raw result rows.
7. Publish the normalized table only after validation and write a checksummed
   completion manifest.

The scDRS path performs these steps:

1. Resolve the selected scDRS v1.0.3 implementation and validate its reported
   version, exact H5AD structure, matrix values/state, annotations, and optional
   covariates in backed chunks.
2. In direct file mode, validate the supplied `.gs`. In pipeline/MAGMA mode,
   read configured gene and Z-statistic columns, validate a one-to-one identifier
   crosswalk, write an auditable gene-statistics table, and call the official
   `scdrs munge-gs` command. Packaged construction requires at least 100 mapped
   input genes and selects up to the top 1,000 with Z-score weights.
3. Validate the generated gene set after intersection with genes surviving the
   configured H5AD filters.
4. Call `scdrs compute-score`; call `perform-downstream` once when any group,
   correlation, or gene analysis is selected.
5. Keep all work in an isolated staging directory. Publish the complete engine
   directory atomically, then checksum inputs and every expected output.

## Outputs

The primary output is the configured `results/<dataset>_magma_celltype.tsv`
table. It contains the dataset identifier, cell type, MAGMA gene count,
coefficient, standardized coefficient, standard error, raw p-value, adjusted
p-values, and significance flags. The native MAGMAcovar result and log are
retained below the configured engine directory. The module also writes its
canonical log, resolved configuration, and completion manifest.

scDRS outputs are retained under `engines/<dataset>_scdrs/`:

- `<trait>.score.gz` and `<trait>.full_score.gz`;
- `<trait>.scdrs_group.<annotation>` for each group annotation;
- `<trait>.scdrs_cell_corr` when continuous correlations are selected;
- `<trait>.scdrs_gene` when gene analysis is selected;
- a PostGWAS QC report containing H5AD, gene-set, covariate, software, and
  scientific-setting provenance, annotation group counts, and the fixed
  upstream scDRS 1.0.3 internal seed;
- for MAGMA construction, the mapped gene-statistics table, generated `.gs`,
  and complete identifier-mapping report.

MAGMA cell typing and scDRS use independent completion manifests, so one
selected method is never mistaken for completion of the other.

## QC and logs

Review the selected primary MAGMA mapping, gene identifier system, gene overlap,
missingness per property, exact MAGMAcovar command, number of tested cell types,
correction family size, adjusted significant counts, executable identity, and
all input/output fingerprints. A completed MAGMAcovar engine subdirectory alone
does not mean the enclosing single-cell result completed; verify the
single-cell completion manifest.

For scDRS, also review the effective post-filter gene-set size and fraction,
matrix-state declaration, cells and genes retained by minimal filtering,
annotation group sizes, covariate overlap/intercept validation, native scDRS
version, internal seed, commands, and the tool-specific completion manifest.

## Interpretation

A positive result indicates that genes with higher expression in a cell type
have stronger gene-level GWAS association after adjustment for average
expression and MAGMA's technical covariates. It does not establish that the cell
type is causal, that the associated genes act only in that cell type, or that
the result is independent of correlated cell-type expression profiles.

A positive scDRS score or group association means that GWAS-prioritized genes
are unusually enriched in the cell's expression profile relative to matched
control genes. It does not identify causal cells or genes and does not provide
the direction of a genetic, expression, abundance, or therapeutic effect.

## Common problems

Runs stop for an absent `Average` column, duplicate or nonnumeric properties,
constant expression profiles, excessive missingness, insufficient gene-ID
overlap, a non-calibrated primary MAGMA mapping, an incomplete MAGMA result, or
an output hypothesis family that differs from the input cell-type columns.

scDRS additionally stops for an unreadable H5AD, duplicate cell/gene IDs,
negative/non-finite expression, fractional values declared as raw counts,
missing annotations, malformed `.gs`, insufficient post-overlap genes, an
oversized effective set, nonnumeric covariates, a crosswalk that is not
one-to-one, a MAGMA identifier-type mismatch, an unsupported scDRS version, or
partial native output. Covariates also fail when they lack the configured
constant column, the constant differs from one, or their cell overlap is not
above the configured threshold.

## Limitations

The MAGMA path covers FUMA workflow step 1 for one aggregated dataset. FUMA's
optional forward-selection and cross-dataset conditional analyses are not yet
implemented. Stratified LDSC, seismic, and scPagwas are not implemented.

The scDRS implementation reads `adata.X` only. It does not yet materialize
`raw.X` or a named layer into a derived H5AD, remove duplicated CELLxGENE
observations, download an atlas/crosswalk, estimate atlas-scale memory, adjust
native p-values, or meta-analyze independent atlases. The upstream scDRS CLI
does not expose a random-seed option. The pinned implementation calls its
scoring method with the fixed internal default seed 0; PostGWAS records both
that seed and that the general PostGWAS seed was not applied. A normalized
matrix declaration is validated for non-negativity and finiteness but cannot
prove which normalization procedure produced the stored values.

Pipeline scDRS uses the resolved general MAGMA configuration and records its
window, ancestry, mapping, and model. The packaged PostGWAS MAGMA window is not
silently relabelled as the original scDRS paper's 10/10-kb profile.

MAGMA text results do not encode enough metadata to prove genome-build,
ancestry, or gene-identifier compatibility in direct mode. Retain the upstream
MAGMA configuration and reference manifest.

### Reference inputs and downloads

PostGWAS deliberately does not select or silently download a biological atlas
or an identifier mapping. Pin the exact files used by the study:

| Input | Appropriate source | What must be checked before use |
|---|---|---|
| H5AD atlas | [CELLxGENE Discover published-data downloads](https://cellxgene.cziscience.com/docs/03__Download%20Published%20Data), an atlas author's repository, or another study archive | RNA assay, tissue/cell subset, primary observations, whether counts are in `X` or `raw.X`, normalization, complete-gene rather than HVG-only content, cell annotations, and consent/use constraints. |
| Entrez-to-symbol crosswalk | A versioned [HGNC archive](https://www.genenames.org/download/archive/), NCBI Gene release, or Ensembl BioMart export | Species/build/release, source and target columns, deprecated IDs, and one-to-one mappings after restricting to the analyzed genes. Preserve the downloaded file and checksum. |
| Native `.gs` | Study-authored gene set or the official [`scdrs munge-gs` format](https://martinjzhang.github.io/scDRS/file_format.html) | Exact `TRAIT`/`GENESET` header, one consistent weighted or unweighted representation per trait, identifier compatibility, and effective overlap after filtering. |
| MAGMA gene result | PostGWAS MAGMA pipeline output, or a separately retained MAGMA run for direct mode | Genome build, ancestry-matched LD, gene-location release, window/model, identifier type, and complete successful provenance. |

Install the pinned implementation and compatible scientific stack with
`pip install -e '.[analysis,single-cell]'` when using a source checkout.

## Scientific references

- [scDRS implementation evidence review](../reference/scdrs-evidence-review.md)
- [Zhang et al. 2022, scDRS](https://doi.org/10.1038/s41588-022-01167-z)
- [FUMA cell-type tutorial and three-step workflow](https://fuma.ctglab.nl/tutorial#cell-type)
- [Watanabe et al. 2019, FUMA cell-type analysis](https://doi.org/10.1038/s41467-019-11181-1)
- [de Leeuw et al. 2015, MAGMA](https://doi.org/10.1371/journal.pcbi.1004219)
- [Official MAGMA documentation](https://cncr.nl/research/magma/)
