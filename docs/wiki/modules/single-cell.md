# Single-Cell Integration

## Purpose

The `single_cell` module provides three complementary integrations of GWAS
evidence with cell-type data: FUMA-compatible MAGMA cell typing, scDRS cell
scoring, and LDSC-SEG cell-type-specific heritability analysis. Select one or
more in execution order with `--tools magma_celltype scdrs ldsc_celltype`.

## What the analysis does

For `magma_celltype`, the pipeline reuses PostGWAS MAGMA to convert
variant-level GWAS evidence into calibrated gene-association results. It then
reuses the MAGMAcovar service to
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

`ldsc_celltype` implements the official LDSC `--h2-cts` workflow. Each `.ldcts`
line defines one annotation to test and one or more control annotations. LDSC
runs a separate stratified LD-score regression for every line, jointly including
the configured baseline model. The reported coefficient is the coefficient of
the first LD-score prefix on that line; PostGWAS preserves its native one-sided
test against a coefficient greater than zero and applies the configured
multiple-testing methods across the complete `.ldcts` hypothesis family.

## When to use it

Use `magma_celltype` when the question is whether cell-type average-expression
properties predict MAGMA gene association. Use `scdrs` when the question is
which individual cells, annotated groups, or cell states preferentially express
a GWAS-derived disease gene set. Use `ldsc_celltype` when the question is whether
SNP heritability is enriched in a predefined cell-type annotation after baseline
and all-genes control adjustment. Use pipeline mode when starting from a
harmonised GWAS-VCF and direct mode with exact method-specific inputs.

## Input requirements

Missing inputs for all selected methods are reported together with their public
CLI options and canonical YAML keys, before reference scans or native commands.
Pipeline-generated gene results, gene sets, and summary statistics are exempt
from this startup requirement; the atlas, identifier crosswalk (when configured),
and LD-score resources are not. A supplied but invalid file still fails the
method's scientific validation and is not relabelled as a missing argument.

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
  Native list-valued options reserve commas. Proportion adjustment additionally
  requires complete group labels and fewer than one group per ten analyzed cells,
  matching the pinned upstream stability check.
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
- The current validator requires identical H5AD and gene-set species values,
  including the same alias spelling. Native cross-species homolog conversion is
  not exposed until PostGWAS can validate the post-conversion effective gene set
  itself; exact spelling also avoids an alias bug in the pinned upstream CLI.

For `ldsc_celltype`:

- Direct mode requires pre-munged LDSC summary statistics containing `SNP`, `Z`,
  and `N`. Pipeline mode obtains the formatter's validated LDSC table and runs
  `munge_sumstats` with the configured INFO/MAF thresholds and an explicit
  HapMap3 merge-alleles file. The pipeline pins the LDSC formatter target to
  rsIDs and rejects a conflicting global or target-specific identifier policy.
- Both modes require a two-field `.ldcts` file. Field one is the unique test
  label. Field two is a comma-delimited list whose first prefix is the tested
  cell-type annotation and whose later prefix or prefixes are controls. The
  official released sets use an all-genes control as the second prefix.
- Relative prefixes inside `.ldcts` are resolved against that file's directory.
  Baseline and weight prefixes supplied by CLI/YAML are resolved against the
  run's working directory. A prefix may contain one `@` chromosome placeholder;
  otherwise LDSC chromosome numbers are appended directly.
- Packaged validation requires chromosomes 1–22. Every baseline and cell-type
  prefix needs nonempty `<prefix><chr>.l2.ldscore.gz` and
  `<prefix><chr>.l2.M_5_50`; the regression-weight prefix needs the LD-score
  files. The code does not guess alternative suffixes or accept a partial panel.
- Summary statistics, baseline scores, weights, cell-type scores, and controls
  must use the same genome build and an ancestry-appropriate LD reference.
  These properties are declared and recorded but cannot be inferred reliably
  from LDSC files, so the user must verify the reference release.

The `.ldcts` file has no header. A minimal relative-prefix example is:

```text
Neuron references/Neuron.,references/AllGenes.
Microglia references/Microglia.,references/AllGenes.
```

The expression covariate matrix and H5AD atlas are processed biological inputs.
The GWAS formatter does not create them, and this module does not perform FASTQ
processing, alignment, study-level cell QC, clustering, or annotation.

## Command

```text
postgwas single_cell --tools magma_celltype --magma-gene-results-file PATH \
  --single-cell-covariates PATH [options]

postgwas single_cell --tools scdrs --scdrs-h5ad-file ATLAS.h5ad \
  --scdrs-gene-set-file TRAIT.gs [options]

postgwas single_cell --tools scdrs --scdrs-gene-set-source magma \
  --scdrs-h5ad-file ATLAS.h5ad \
  --scdrs-magma-gene-results-file STUDY.genes.out \
  --scdrs-gene-id-map entrez_to_symbol.tsv [options]

postgwas single_cell --tools ldsc_celltype \
  --ldsc-celltype-sumstats-file STUDY.sumstats.gz \
  --ldsc-celltype-ldcts-file BRAIN.ldcts \
  --ldsc-celltype-baseline-prefix reference/baselineLD. \
  --ldsc-celltype-weights-prefix reference/weights. [options]
```

## Direct mode

### MAGMA cell typing

```console
postgwas single_cell \
  --tools magma_celltype \
  --magma-gene-results-file STUDY.genes.raw \
  --single-cell-covariates atlas_celltype_average.tsv \
  --dataset-id STUDY \
  --output-directory results
```

MAGMA is resolved from `resources.executables.magma` (default `magma`) and
validated on `PATH` before analysis.

### scDRS from an existing gene set

This example assumes human gene symbols in the `.gs` and H5AD, raw counts in
`adata.X`, and a `cell_type` annotation in `adata.obs`. Declare the actual matrix
state with `--scdrs-matrix-state`; do not label log-normalized values as counts.

```console
postgwas single_cell \
  --tools scdrs \
  --scdrs-h5ad-file brain_atlas.h5ad \
  --scdrs-gene-set-file STUDY.gs \
  --scdrs-matrix-state raw_counts \
  --scdrs-group-analysis cell_type \
  --dataset-id STUDY \
  --output-directory results
```

### scDRS from existing MAGMA gene statistics

This direct mode constructs the disease gene set from a completed, headered
`.genes.out`, not the `.genes.raw` used by MAGMA cell typing. The example assumes
Entrez MAGMA identifiers, a one-to-one tab-delimited crosswalk with `ENTREZID`
and `SYMBOL` headers, human gene symbols in `adata.var_names`, and raw counts
in `adata.X`. Keep the crosswalk release and MAGMA gene universe consistent.

```console
postgwas single_cell \
  --tools scdrs \
  --scdrs-gene-set-source magma \
  --scdrs-magma-gene-results-file STUDY.genes.out \
  --scdrs-h5ad-file brain_atlas.h5ad \
  --scdrs-gene-id-map entrez_to_symbol.tsv \
  --scdrs-matrix-state raw_counts \
  --scdrs-group-analysis cell_type \
  --dataset-id STUDY \
  --output-directory results
```

### LDSC cell typing from munged statistics

Use a matched GRCh37/EUR reference release for this example. The direct input
is already munged; a formatter table is not interchangeable with it.

```console
postgwas single_cell \
  --tools ldsc_celltype \
  --ldsc-celltype-sumstats-file STUDY.sumstats.gz \
  --ldsc-celltype-ldcts-file reference/brain.ldcts \
  --ldsc-celltype-baseline-prefix reference/baselineLD. \
  --ldsc-celltype-weights-prefix reference/weights. \
  --ldsc-celltype-genome-build GRCh37 \
  --ldsc-celltype-population EUR \
  --dataset-id STUDY \
  --output-directory results
```

## Pipeline mode

### MAGMA cell typing

This example uses the packaged GRCh37/EUR MAGMA settings and NCBI/Entrez gene
locations. The expression matrix must use matching Entrez IDs. For an Ensembl
atlas, supply compatible Ensembl gene locations and their
[primary-ID and source declarations](magma.md#gene-identifiers-and-source-declarations);
changing an identifier label does not convert either file.

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

### scDRS

For pipeline scDRS, the same stages are used, but the final stage consumes
MAGMA's headered gene output and creates a weighted disease gene set first:

```console
postgwas pipeline \
  --modules single_cell \
  --tools scdrs \
  --vcf study.gwas.vcf.gz \
  --scdrs-h5ad-file brain_atlas.h5ad \
  --scdrs-gene-id-map entrez_to_symbol.tsv \
  --scdrs-matrix-state raw_counts \
  --scdrs-group-analysis cell_type \
  --magma-ld-reference reference/g1000_eur \
  --gene-location-file reference/NCBI37.3.gene.loc \
  --dataset-id STUDY \
  --output-directory results
```

The packaged crosswalk schema expects `ENTREZID` and `SYMBOL`; the H5AD example
assumes raw counts in `adata.X` and gene symbols in `adata.var_names`. The
namespace declarations can be overridden with `--scdrs-source-gene-id-type` and
`--scdrs-target-gene-id-type`, but mapping-column names and mapping policy are
configured in YAML under `modules.single_cell.scdrs.magma_gene_set`. Keep these
settings consistent with the actual crosswalk. Set `mapping_mode: exact` only
when MAGMA and `adata.var_names` genuinely use the same identifiers.

### LDSC cell typing

For LDSC cell typing, pipeline mode schedules `formatter` and `single_cell`
without MAGMA. It creates the LDSC input from the harmonised GWAS-VCF, mungs it,
then runs the cell-type regressions:

```console
postgwas pipeline \
  --modules single_cell \
  --tools ldsc_celltype \
  --vcf study.gwas.vcf.gz \
  --ldsc-celltype-ldcts-file reference/brain.ldcts \
  --ldsc-celltype-baseline-prefix reference/baselineLD. \
  --ldsc-celltype-weights-prefix reference/weights. \
  --ldsc-celltype-merge-alleles-file reference/w_hm3.snplist \
  --ldsc-celltype-genome-build GRCh37 \
  --ldsc-celltype-population EUR \
  --ldsc /path/to/ldsc.py \
  --munge-sumstats /path/to/munge_sumstats.py \
  --dataset-id STUDY \
  --output-directory results
```

When LDSC cell typing is selected together with a MAGMA-based method, the
planner schedules the union once: `formatter`, `magma`, and `single_cell`.

## Parameters

`modules.single_cell.tools` is the canonical ordered tool selection and is
overridden by `--tools`. The accepted values are `magma_celltype`, `scdrs`, and
`ldsc_celltype`.
Each tool validates only its own inputs. MAGMA cell typing retains its normalized
result schema and multiple-testing family. scDRS has separate H5AD, gene-set,
mapping, filtering, control-set, downstream-analysis, version, and native-output
settings under `modules.single_cell.scdrs`.

LDSC cell-type inputs, reference-file syntax, output columns, native suffixes,
declared build/population, and version probe are under
`modules.single_cell.ldsc_celltype`. The reused formatter-to-LDSC munging
thresholds remain under `modules.ldsc`; resolved metadata and exported pipeline
configuration retain both configuration sections.

`magma_celltype.magmacovar_use_case` selects a named definition from
`modules.magmacovar.model_use_cases`. PostGWAS verifies that the selected use
case resolves to the FUMA base model and positive one-sided direction. The
remaining MAGMAcovar missingness, overlap, executable, staging, and validation
settings continue to be owned and logged by the reused MAGMAcovar service.

## Processing steps

The `magma_celltype` path performs these steps:

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

The LDSC cell-type path performs these steps:

1. Resolve the LDSC executable and version; parse `.ldcts` exactly; reject blank,
   malformed, duplicate, or repeated tested-annotation entries.
2. Resolve every prefix and require the complete configured autosomal LD-score
   and regression-SNP-count inventory before starting analysis.
3. In direct mode, validate the pre-munged `SNP/Z/N` table. In pipeline mode,
   reuse the formatter's LDSC table and existing LDSC munging command/settings,
   then validate the generated `.sumstats.gz`.
4. Write a normalized absolute-prefix `.ldcts`, then invoke native LDSC with
   `--h2-cts`, baseline `--ref-ld-chr`, regression `--w-ld-chr`, and
   `--ref-ld-chr-cts`.
5. Require the native result labels to match the complete input family, validate
   finite coefficients, positive standard errors, and P values in `[0,1]`, then
   apply the configured correction methods without replacing the native P value.
6. Publish the engine and normalized result from an isolated partial directory,
   write QC/reference provenance, and checksum the inputs and outputs for resume.

## Outputs

For `magma_celltype`, the primary output is the configured
`results/<dataset>_magma_celltype.tsv` table. It contains the dataset identifier,
cell type, MAGMA gene count, coefficient, standardized coefficient, standard
error, raw p-value, adjusted
p-values, and significance flags. The native MAGMAcovar result and log are
retained below the configured engine directory. The module also writes its
canonical log, resolved configuration, and completion manifest.

scDRS outputs are retained under `engines/<dataset>_scdrs/`:

- `<trait>.score.gz` and `<trait>.full_score.gz`;
- `<trait>.scdrs_group.<annotation>` for each group annotation;
- `<trait>.scdrs_cell_corr` when continuous correlations are selected;
- `<trait>.scdrs_gene` when gene analysis is selected;
- a PostGWAS QC report containing H5AD, gene-set, covariate, software, and
  scientific-setting provenance, input and post-filter analysis group counts,
  and the fixed upstream scDRS 1.0.3 internal seed;
- for MAGMA construction, the mapped gene-statistics table, generated `.gs`,
  and identifier-mapping summary; the checksummed source crosswalk remains the
  row-level mapping provenance.

All three methods use independent, tool-scoped completion manifests. Changing
an input owned by one method does not invalidate an otherwise identical result
from another, and one selected method is never mistaken for completion of a
different method.

LDSC cell typing writes:

- `results/<dataset>_ldsc_celltype.tsv`, using the shared normalized cell-type
  schema. `beta` is the native LDSC annotation coefficient; gene-count and
  standardized-beta fields are `NA` because LDSC does not report them;
- `engines/<dataset>_ldsc_celltype/<dataset>.cell_type_results.txt` and the
  native LDSC `.log`;
- the absolute-prefix `.ldcts` actually passed to LDSC;
- in pipeline/formatter mode, the munged `.sumstats.gz` and munging log;
- a QC YAML containing software version, cell types, declared build/ancestry,
  scientific interpretation, exact reference inventory with sizes and
  modification times, and summary metrics;
- a tool-scoped completion manifest. Reference panels are represented in its
  configuration digest by the same inventory instead of hashing every large LD
  file; summary statistics, `.ldcts`, merge list, and produced files are hashed.

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

For LDSC cell typing, review the declared genome build and population, all 22
validated chromosomes, baseline/weight/cell-type/control prefixes, normalized
`.ldcts`, HapMap3 merge list, munging exclusions in the native log, LDSC version,
coefficient/standard error, one-sided P value, complete correction family, and
tool-specific completion manifest.

## Interpretation

For `magma_celltype`, a positive result indicates that genes with higher
expression in a cell type have stronger gene-level GWAS association after adjustment for average
expression and MAGMA's technical covariates. It does not establish that the cell
type is causal, that the associated genes act only in that cell type, or that
the result is independent of correlated cell-type expression profiles.

A positive scDRS score or group association means that GWAS-prioritized genes
are unusually enriched in the cell's expression profile relative to matched
control genes. It does not identify causal cells or genes and does not provide
the direction of a genetic, expression, abundance, or therapeutic effect.

A positive LDSC cell-type result means the tested annotation has a positive
additional per-SNP heritability contribution conditional on the baseline model
and the other annotation prefixes on that `.ldcts` line. It is not total cell-
type heritability, a causal cell-type probability, a locus/gene list, or evidence
that the result is independent of correlated annotations. The native P value is
one-sided; use a configured adjusted P value for the tested family.

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

LDSC cell typing additionally stops for malformed `.ldcts`, fewer than two
prefixes per line, repeated tested annotations, missing baseline/weight/control
prefixes, any missing chromosome LD score or `.l2.M_5_50`, missing `SNP/Z/N`,
unparseable LDSC version, a native result family that differs from `.ldcts`,
non-finite coefficients, nonpositive standard errors, partial native output, or
provenance mismatch during resume.

## Limitations

The MAGMA path covers FUMA workflow step 1 for one aggregated dataset. FUMA's
optional forward-selection and cross-dataset conditional analyses are not yet
implemented. LDSC `--h2-cts` is implemented, but general-purpose annotation
construction/LD-score estimation and arbitrary partitioned-heritability models
are not part of this module. CELLECT-LDSC, SEISMIC, and scPagwas are not
implemented.

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

LDSC LD-score files likewise do not carry enough standardized metadata to prove
their genome build, population, annotation derivation, or compatibility with the
GWAS. PostGWAS verifies the declared, complete file contract and records it; the
researcher must pin and review the release. The method consumes precomputed
cell-type LD scores and does not derive annotations or LD scores from H5AD.

### Reference inputs and downloads

PostGWAS deliberately does not select or silently download a biological atlas
or an identifier mapping. Pin the exact files used by the study:

| Input | Appropriate source | What must be checked before use |
|---|---|---|
| H5AD atlas | [CELLxGENE Discover published-data downloads](https://cellxgene.cziscience.com/docs/03__Download%20Published%20Data), an atlas author's repository, or another study archive | RNA assay, tissue/cell subset, primary observations, whether counts are in `X` or `raw.X`, normalization, complete-gene rather than HVG-only content, cell annotations, and consent/use constraints. |
| Entrez-to-symbol crosswalk | A versioned [HGNC archive](https://www.genenames.org/download/archive/), NCBI Gene release, or Ensembl BioMart export | Species/build/release, source and target columns, deprecated IDs, and one-to-one mappings after restricting to the analyzed genes. Preserve the downloaded file and checksum. |
| Native `.gs` | Study-authored gene set or the official [`scdrs munge-gs` format](https://martinjzhang.github.io/scDRS/file_format.html) | Exact `TRAIT`/`GENESET` header, one consistent weighted or unweighted representation per trait, identifier compatibility, and effective overlap after filtering. |
| MAGMA gene result | PostGWAS MAGMA pipeline output, or a separately retained MAGMA run for direct mode | Genome build, ancestry-matched LD, gene-location release, window/model, identifier type, and complete successful provenance. |
| Cell-type `.ldcts` plus LD scores | A versioned [official LDSC-SEG release](https://github.com/bulik/ldsc/wiki/Cell-type-specific-analyses), or LD scores created from a documented annotation and ancestry-matched reference panel | The first prefix is the tested annotation, later prefixes are controls, chromosomes 1–22 are complete, `.l2.M_5_50` files match the LD scores, and annotation/build/population provenance is retained. |
| Baseline LD scores, weights, and HapMap3 list | The matching release linked by the [official LDSC cell-type tutorial](https://github.com/bulik/ldsc/wiki/Cell-type-specific-analyses) | All resources belong to the same reference population/build/release and use the SNP universe expected by the selected LDSC implementation. Do not mix baseline, weights, and cell-type releases. |

For the complete software setup, follow the
[installation guide](../getting-started/installation.md). In a source checkout,
`pip install -e '.[analysis,single-cell]'` installs the declared Python extras;
it does not by itself install the MAGMA/LDSC executables or download the atlases,
crosswalks, and LD-score reference panels described above.

### Implementation architecture

The module uses a registry-backed method hierarchy so each integration method
owns its command-line options, preflight, execution, validation, and result
publication:

```text
single_cell/
├── cli.py
├── service.py
└── methods/
    ├── base.py
    ├── registry.py
    ├── magma_celltype/
    │   ├── cli.py
    │   ├── analysis.py
    │   └── method.py
    ├── scdrs/
    │   ├── cli.py
    │   ├── method.py
    │   └── runner.py
    └── ldsc_celltype/
        ├── cli.py
        ├── analysis.py
        ├── method.py
        └── runner.py
```

The shared service resolves configuration once, executes registered methods in
the order given by `modules.single_cell.tools`, and combines non-colliding
artifacts. MAGMA cell typing delegates native execution to the existing
`magmacovar` module. LDSC cell typing reuses the existing formatter's LDSC
contract and the existing LDSC command builders and munging policy. New methods
must implement the common method contract and be added to the registry and
schema; all method defaults remain in the canonical `single_cell.yaml`.

## Scientific references

- [scDRS implementation evidence review](../reference/scdrs-evidence-review.md)
- [Zhang et al. 2022, scDRS](https://doi.org/10.1038/s41588-022-01167-z)
- [FUMA cell-type tutorial and three-step workflow](https://fuma.ctglab.nl/tutorial#cell-type)
- [Watanabe et al. 2019, FUMA cell-type analysis](https://doi.org/10.1038/s41467-019-11181-1)
- [de Leeuw et al. 2015, MAGMA](https://doi.org/10.1371/journal.pcbi.1004219)
- [Official MAGMA documentation](https://cncr.nl/research/magma/)
- [Official LDSC cell-type-specific analysis workflow](https://github.com/bulik/ldsc/wiki/Cell-type-specific-analyses)
- [Finucane et al. 2018, LDSC-SEG](https://doi.org/10.1038/s41588-018-0081-4)
