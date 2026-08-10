# scDRS Evidence Review

## Purpose and review scope

This page records the scientific evidence used to implement and review scDRS in
the PostGWAS `single_cell` module. It separates three kinds of
information:

- **Normative evidence**: assumptions and behavior stated by the scDRS authors,
  current official documentation, and upstream source code.
- **Application evidence**: choices reported by studies that applied scDRS to
  real GWAS and single-cell data.
- **PostGWAS policy**: proposed validation and provenance requirements. These
  are design decisions, not behavior implemented by the current release.

Every paper listed in the evidence tables below was reviewed for the method
details reported here. This is an implementation-focused review, not a claim to
have reviewed every publication that merely cites or mentions scDRS. Papers
were included when they described a distinct scDRS input, preprocessing,
parameter, validation, scaling, or interpretation decision. Reviews and
non-RNA adaptations were excluded from the implementation evidence.

Evidence and upstream software status were checked on **2026-08-09**. PostGWAS
supports exact-H5AD/exact-`.gs` direct mode and MAGMA-to-`.gs` pipeline mode
with an explicit one-to-one identifier crosswalk.

## Normative scDRS contract

The primary method paper defines scDRS as a GWAS-to-single-cell approach. It
first constructs a putative disease-gene set, then compares its aggregate
expression in each cell with expression of matched control gene sets. Matching
accounts for gene-set size and the expression mean and variance of the disease
genes. The method returns individual-cell scores and P values; downstream tests
summarize association and heterogeneity for annotated cell groups.

The following table is the implementation contract supported by the
[primary publication](https://doi.org/10.1038/s41588-022-01167-z),
[official CLI reference](https://martinjzhang.github.io/scDRS/reference_cli.html),
[file-format reference](https://martinjzhang.github.io/scDRS/file_format.html),
and [upstream repository](https://github.com/martinjzhang/scDRS).

| Topic | Upstream requirement or recommendation | Consequence for PostGWAS |
|---|---|---|
| Expression matrix | `adata.X` is assumed to be size-factor normalized, for example to 10,000 counts per cell, and transformed with `log(x + 1)`. | Matrix source and preprocessing state must be declared; PostGWAS must not guess them. |
| Raw-count mode | `--flag-raw-count True` applies size-factor normalization and `log1p`; `False` assumes this has already been done. | Double normalization must fail validation. Raw counts stored outside `X` require a prepared copy. |
| Minimal filtering | CLI defaults are at least 250 detected genes per cell and expression in at least 50 cells per gene. | These are scDRS minimal filters, not a replacement for study-level cell QC. Effective values must be configured and logged. |
| Numerical values | The loader rejects `NaN` and negative expression values. | Scaled or residualized matrices containing negative values are not valid scDRS expression input. |
| Disease-gene set | The method assumptions require a moderate post-overlap set, described as more than 50 genes and less than 20% of the expression-gene universe. | Validation must occur after identifier conversion and expression overlap, not only on the original `.gs` file. |
| Default GWAS gene set | The publication used the top 1,000 MAGMA genes, weighted by MAGMA Z score. Current `munge-gs` defaults permit 100 to 1,000 genes and use Z-score weights. | Top-1,000 weighted MAGMA genes should be a reproducible profile, not an unchangeable algorithm invariant. |
| MAGMA mapping | The publication used a 10-kb window around gene bodies and a population-matched LD reference. | A MAGMA result may be reused only when its build, LD population, gene identifiers, gene window, and model match the selected scDRS profile. |
| GWAS power | The authors recommend an LDSC heritability Z score above 5 or GWAS sample size above 100,000 to obtain a reasonable number of discoveries. | This is a power recommendation rather than a validity threshold. Report it prominently and warn when neither condition is met. |
| Cell diversity | The authors recommend a diverse set of potentially relevant cells; they state that fewer cells should reduce power rather than create false positives. | Dataset relevance and composition must be reported. No universal tissue or disease-status filter can be inferred. |
| Control gene sets | Default `n_ctrl` is 1,000. Runtime and memory scale linearly with cells and control sets. | `n_ctrl`, random seed, expected P-value resolution, runtime, and memory must be resolved before execution. |
| Cell-group proportions | `adj_prop` is recommended only for extremely unbalanced datasets. | Do not enable proportion adjustment automatically; validate the named annotation and record why it was selected. |
| Covariates | Covariates are numerical, indexed by `adata.obs_names`, include a constant, and must have values for at least 75% of cells. | Covariate inclusion is a study-level decision. Missingness, rank, alignment, and transformations require preflight validation. |
| Group analysis | The group-association statistic uses the upper 5% quantile; heterogeneity uses Geary's C. | The cell-group annotation must be present and cell counts must be reported before testing. |
| P values | Cell and group result files report raw P values; the number of controls limits Monte Carlo resolution. | PostGWAS must preserve raw results and define an explicit multiple-testing family rather than relabel raw values as adjusted. |
| Scale | The paper reports about 3 hours and 60 GB for one million cells; a later brain application reported approximately 40 hours and 600 GB for 3.3 million nuclei. | Backed inspection and a pre-run resource estimate are required for atlas-scale data. |

The upstream README currently calls v1.0.3b the latest stable version and v1.0.4
a development version. The online API documentation is labelled v1.0.4.
PostGWAS therefore must pin and record a supported executable version instead
of assuming that the documentation version and installed package are the same.

## CELLxGENE input contract

CELLxGENE Discover standardizes AnnData structure and a minimum metadata
schema; it does not apply one uniform biological QC and normalization workflow
to every submitted study. The
[current CELLxGENE schema](https://chanzuckerberg.github.io/single-cell-curation/latest-schema.html)
requires raw RNA measurements when a standard raw matrix exists and strongly
recommends a normalized matrix. For scRNA-seq, raw counts are normally in
`adata.raw.X` when normalized data are present in `adata.X`; otherwise raw counts
may be in `adata.X`. The normalized method is not standardized across studies.

CELLxGENE gene indices are stable reference identifiers, normally Ensembl gene
IDs without version suffixes, while `var["feature_name"]` supplies the
human-readable gene name. Cell, tissue, disease, organism, assay, donor, and
`is_primary_data` metadata follow controlled schema fields, but labels may be
`unknown`, and a pre-analysis dataset may not yet provide the annotations
needed for scDRS group analysis. Published downloads are documented as
schema-compliant `.h5ad` files in the
[CELLxGENE download guide](https://cellxgene.cziscience.com/docs/03__Download%20Published%20Data).

Consequently, PostGWAS must not pass a CELLxGENE download directly to scDRS
without resolving all of the following:

1. RNA assay and biological subset to analyze.
2. `X`, `raw.X`, or a named layer as the expression source.
3. Raw-count versus size-factor-normalized and `log1p` state.
4. Primary versus duplicated observations.
5. Gene identifier source and an unambiguous mapping to the `.gs` identifiers.
6. Required cell-group annotations, missing labels, and group sizes.
7. Donor, batch, disease-status, and other candidate covariates.

The downloaded artifact must remain unchanged. Any matrix conversion must be
written as a new, provenance-linked `.h5ad` file.

## Published implementation evidence

### Primary method and software

| Study or source | Data and implementation | Key finding for implementation |
|---|---|---|
| [Zhang et al. 2022](https://doi.org/10.1038/s41588-022-01167-z) | 74 traits and 16 scRNA-seq datasets; top 1,000 MAGMA genes weighted by Z score; 1,000 matched controls; size-factor-normalized and log-transformed expression. | Defines the reference workflow, gene-set-size assumptions, GWAS power guidance, diverse-cell recommendation, optional covariate regression, and proportion-adjustment caution. |
| [scDRS CLI](https://martinjzhang.github.io/scDRS/reference_cli.html) | `munge-gs`, `compute-score`, and `perform-downstream`; filtering defaults 250 genes per cell and 50 cells per gene; raw-count normalization flag; default 1,000 controls. | These options must be mapped to one schema-validated PostGWAS configuration and recorded exactly. |
| [scDRS loader](https://martinjzhang.github.io/scDRS/reference/scdrs.util.load_h5ad.html) | Reads `adata.X`, optionally filters and normalizes, and rejects `NaN` or negative expression. | PostGWAS must explicitly materialize the selected CELLxGENE matrix in `X`; native scDRS does not select `raw.X` on the user's behalf. |

### Applications and benchmarks

| Study | Reported choices | Evidence contributed to the PostGWAS design |
|---|---|---|
| [Mosquera et al. 2023, human atherosclerosis](https://doi.org/10.1016/j.celrep.2023.113380) | scDRS v1.0.2; top 1,000 MAGMA genes with Z weights; raw counts exported to H5AD and processed with Scanpy; sex and disease status covariates; 250 control sets; analysis across the full diverse meta-atlas. | Supports raw-count reconstruction, explicit covariates, and diverse cells, but demonstrates that applications sometimes reduce controls below the reference 1,000. |
| [Siewert et al. 2024, orofacial clefting](https://doi.org/10.1038/s41598-024-77724-9) | Removed scaled data before H5AD conversion; retained counts and normalized expression; tested an 87-gene GWAS/TAD set with MAGMA weights and an unweighted sensitivity analysis; correlated score with detected molecules. | Shows a valid custom gene set above the method's lower-size guidance and provides a useful technical-bias diagnostic. It argues against using scaled or integrated values as expression input. |
| [Townsend et al. 2024, inflammatory-disease benchmark](https://doi.org/10.3389/fimmu.2024.1454263) | scDRS v1.0.2 defaults; MAGMA P values passed to `munge-gs`; covariates included nUMI, detected genes, sex, and study-specific age, duration, location, mitochondrial fraction, or smoking; evaluated 100, 300, 500, 1,000, 1,500, and 2,000 top genes. | Demonstrates that covariates must be dataset-specific and that gene-set-size sensitivity should be available for important results. |
| [Lee et al. 2024 preprint, human microglia and AD](https://doi.org/10.1101/2023.10.25.23297558) | MAGMA GRCh38, 1000 Genomes EUR, 35-kb upstream/10-kb downstream; top 200, 500, and 1,000 gene sensitivity; scDRS bins 20 by 20; 100 controls; separate cohorts followed by cell-count-weighted Stouffer meta-analysis. | Confirms that published MAGMA windows and control counts vary. Independent-cohort concordance is stronger evidence than selecting a parameter in one dataset. Preprint evidence is not normative. |
| [Mathys et al. 2024, multiregion AD brain](https://doi.org/10.1038/s41586-024-07606-7) | Applied MAGMA-based scDRS to approximately 1.3 million nuclei and summarized the fraction of FDR-significant cells by cell type, subtype, and region. | Supports multi-level annotations and region-aware reporting, but the main paper does not disclose enough scDRS preprocessing detail to define defaults. |
| [Zhang et al. 2025, intervertebral-disc degeneration](https://doi.org/10.1002/advs.202500505) | scDRS v1.0.1; MAGMA v1.10; hg19; East Asian 1000 Genomes reference; 10-kb extensions; top 1,000 Z-weighted genes. | Reinforces ancestry-matched LD and shows the reference top-1,000 workflow outside European GWAS. |
| [Zhou and Xue 2025, respiratory diseases](https://doi.org/10.3390/biology14121765) | scDRS v1.0.2 and MAGMA v1.10; top 1,000 disease and control genes; other parameters described as defaults; group association and Geary's-C heterogeneity with FDR below 0.05. | Supports retaining both association and heterogeneity outputs, while illustrating that “defaults” alone is insufficient provenance for PostGWAS. |
| [Mandelia et al. 2025, renal-cell carcinoma](https://doi.org/10.1038/s42003-025-09297-w) | scDRS v1.0.1; MAGMA v1.10; top 1,000 genes; approximately 500,000 cells from normal Tabula Sapiens tissues. | Demonstrates that disease-case cells are not required: scDRS maps inherited risk to relevant expression contexts in a normal reference atlas. |
| [Lai et al. 2025, seismic benchmark](https://doi.org/10.1038/s41467-025-63753-z) | scDRS defaults of 1,000 genes and 1,000 controls; removed groups with fewer than 20 cells for method comparison; converted Entrez identifiers to symbols and averaged Z scores for multiple mappings; used analytical P values derived from Monte Carlo Z scores followed by FDR because empirical resolution was limited. | Identifies gene-mapping ambiguity, rare-group handling, and Monte Carlo resolution as explicit design decisions. Averaging ambiguous mappings is an application choice, not a safe silent default. |
| [Li et al. 2025 preprint, method benchmark](https://doi.org/10.1101/2025.05.24.25328275) | Compared scDRS using MAGMA and mBAT-combo gene statistics; used 10-kb gene extensions and 1000 Genomes EUR; excluded groups with fewer than 20 cells for consistent comparison. | Supports a future alternative gene-statistics input, but the initial PostGWAS implementation should reproduce MAGMA first. Preprint conclusions require caution. |
| [Yao et al. 2025, psychiatric traits and brain regions](https://doi.org/10.1038/s41467-024-55611-1) | Reported that one trait on 3.3 million nuclei required roughly 40 compute hours and 600 GB, making broad scDRS analysis infeasible without down-sampling. | Requires resource estimation and rejects an assumption that every whole-atlas H5AD is practical to score directly. Down-sampling would change the scientific input and must never be automatic. |
| [Li et al. 2026, vascular smooth-muscle states and CAD](https://doi.org/10.1038/s41467-026-70530-z) | Used MAGMA Z weights with scDRS and added an expression-direction analysis because scDRS enrichment itself does not provide effect direction. | Establishes an interpretation boundary: a positive scDRS association is not evidence that increasing expression, cell abundance, or cell activity increases disease risk. |

## Agreements across the evidence

The strongest repeated practices are:

1. Derive gene-level GWAS evidence with MAGMA using a build- and
   ancestry-compatible LD reference.
2. Use Z-score-weighted disease genes; top 1,000 is the most common reference
   profile, not the only scientifically used profile.
3. Score non-negative, unscaled expression that represents all retained genes,
   not a highly-variable-gene-only matrix.
4. Preserve cell-level scores and control scores so association,
   heterogeneity, continuous-variable correlation, and gene correlation can be
   reproduced.
5. Treat sequencing depth, detected genes, batch, donor, and biological
   variables as study-specific covariate candidates rather than a universal
   hardcoded set.
6. Correct downstream hypotheses explicitly and report group size and Monte
   Carlo resolution.
7. Validate important findings across gene-set sizes, datasets, cohorts, or
   related methods when feasible.

## Material disagreements and risks

| Decision | Observed variation | Risk if hidden |
|---|---|---|
| MAGMA gene window | Original scDRS and several applications use 10 kb on both sides; other studies and the current PostGWAS MAGMA default use 35 kb upstream and 10 kb downstream. | The selected genes and weights change. A pre-existing `.genes.raw` file is not interchangeable across windows. |
| Gene-set selection | Top 1,000 is common; studies also use significance thresholds, 87 curated genes, or sensitivity ranges from 100 to 2,000. | Results may be driven by an arbitrary rank cutoff or violate the moderate-size assumption after overlap. |
| Number of controls | 100, 250, and 1,000 are reported. | Group-level empirical P-value resolution is approximately `1 / (n_ctrl + 1)` and may be inadequate after multiple testing. |
| Matrix preparation | Raw counts with scDRS normalization, author-normalized values, integrated objects with scaled values removed, and even imputed values appear in applications. | Double normalization, negative values, denoising artifacts, and loss of the mean-variance relationship can alter calibration. |
| Gene identifiers | MAGMA commonly emits Entrez identifiers; CELLxGENE uses Ensembl indices; scDRS gene sets commonly use symbols. | Silent one-to-many conversion changes the gene set and expression matrix. |
| Covariates | None, technical metrics, demographics, disease status, and sample location have all been used. | Under-adjustment leaves technical bias; over-adjustment can remove the biological signal of interest. |
| Cell-group size | Native scDRS and comparative studies use different minimums; several benchmarks standardized on 20 cells. | Very small groups yield unstable association and heterogeneity summaries. |
| Significance | Native Monte Carlo P values, Monte Carlo Z scores converted to analytical P values, and differing FDR families are reported. | Nominal and adjusted evidence can be confused, especially when control-set resolution is coarse. |
| Atlas scale | Published runs range from tens of thousands to millions of cells. | Unbounded execution can exhaust memory or lead users to down-sample without recording the changed estimand. |

## Proposed PostGWAS scientific policy

The following policy is the recommended starting point for implementation. It
must remain visibly marked as proposed until code, schema, tests, and user
documentation implement the same behavior.

### 1. Reuse MAGMA conditionally

PostGWAS should reuse the existing MAGMA service and a completed MAGMA artifact
only when the following provenance matches the requested scDRS analysis:

- genome build;
- LD-reference population and reference version;
- gene-location reference and identifier type;
- upstream and downstream gene windows;
- gene model;
- study sample-size field and variant universe;
- successful completion manifest.

The reproducibility profile matching the original scDRS publication should use
10-kb upstream and downstream windows, `snp-wise=mean`, top 1,000 genes,
Z-score weights, and 1,000 controls. The current general PostGWAS MAGMA default
of 35/10 kb must not be relabelled as the original scDRS workflow. If the
requested profile differs from the existing MAGMA output, the same MAGMA
service should run a separate, provenance-distinct mapping rather than reuse an
incompatible file.

### 2. Prepare H5AD without changing the source

Direct mode should accept an exact H5AD plus either an exact `.gs` file or an
exact gene-statistics table. Pipeline mode should obtain the gene statistics
from the compatible MAGMA stage. In both modes, configuration must declare:

- matrix source: `X`, `raw.X`, or an exact layer;
- matrix state: raw counts or size-factor-normalized `log1p` values;
- species and assay;
- gene-identifier source;
- group and continuous annotations;
- cell-selection predicates;
- covariate columns or exact covariate file.

For a CELLxGENE RNA download, the preferred reproducible path is to copy raw
counts from `raw.X` into `X` in a derived H5AD and allow scDRS to normalize it.
If `raw.X` is absent, PostGWAS should accept `X` only when the user explicitly
declares its preprocessing state and the validation report supports that
declaration. PostGWAS should not use scaled, centered, batch-corrected,
HVG-only, or imputed expression as the standard profile.

### 3. Make identifier conversion auditable

The conversion from MAGMA Entrez identifiers and CELLxGENE Ensembl identifiers
to the identifiers used by scDRS must use a pinned, build- and species-specific
crosswalk. The report must contain input, mapped, unmapped, ambiguous,
deprecated, duplicate, and post-overlap counts. Ambiguous mappings should fail
by default; neither selecting the first identifier nor appending suffixes makes
a scientifically equivalent gene.

### 4. Validate the effective gene set

Validation must be applied after identifier conversion and intersection with
the expression matrix. The standard profile should fail when the effective set
does not satisfy the original method's moderate-size assumption. A future
non-standard override, if provided, must be explicit and mark the result as
departing from the reference method. Gene-set sizes used for sensitivity
analysis are separate hypothesis-generating runs and must not overwrite the
primary result.

### 5. Keep covariate decisions explicit

PostGWAS should propose candidate columns from metadata but never silently add
them to the model. Technical depth and detected-gene metrics are common
candidates. Disease status, age, sex, donor, region, and batch can be either
confounders or part of the biological question. The resolved configuration and
log must state why each covariate was included, its encoding, missingness, and
whether it is constant or collinear.

### 6. Protect inference and interpretation

PostGWAS should retain native raw scores, normalized scores, control scores,
Monte Carlo P values, and Monte Carlo Z scores. It should calculate adjusted
values only in a declared family such as all cell groups for one trait and one
annotation. Before running, it should report whether `n_ctrl` permits the
requested empirical significance resolution. An analytical P value derived
from a Monte Carlo Z score, if supported, must be a separately named optional
result and never replace the native empirical P value.

Results must be described as enrichment of GWAS-prioritized gene expression.
They do not establish causal cell types, causal genes, effect direction,
disease-state differential expression, or therapeutic direction.

### 7. Require reproducibility artifacts

At minimum, a completed run should publish:

- the prepared H5AD fingerprint and preparation report;
- H5AD structural, matrix, annotation, and group-count QC;
- MAGMA provenance and the exact gene-statistics table;
- identifier-conversion and gene-overlap reports;
- generated `.gs` file and selection report;
- scDRS command, version, seed, resource usage, and native outputs;
- adjusted downstream result tables with declared testing families;
- resolved configuration and a checksummed completion manifest.

## Remaining implementation questions

The initial implementation resolved the executable version as pinned stable
v1.0.3b, rejects ambiguous mappings, and records rather than relabels the
general MAGMA profile. The following extensions remain open:

1. Should an explicitly named original 10/10-kb MAGMA profile be added alongside
   the recorded general PostGWAS MAGMA profile?
2. Which downloadable, checksum-pinned crosswalk bundles should support Entrez,
   Ensembl, and HGNC symbols for
   GRCh37 and GRCh38?
3. What minimum cell-group size is required for association and Geary's-C
   heterogeneity, separately?
4. Is native Monte Carlo inference the only primary result, or will an optional
   analytical-P profile be supported for large testing families?
5. How should independent H5AD datasets be handled: separate runs only, or a
   later, separately specified meta-analysis stage?
6. Which resource estimator and upper limits are required before accepting a
   whole CELLxGENE atlas?

## Reviewed references

- Zhang MJ, Hou K, Dey KK, et al. 2022.
  [Polygenic enrichment distinguishes disease associations of individual cells in single-cell RNA-seq data](https://doi.org/10.1038/s41588-022-01167-z).
- de Leeuw CA, Mooij JM, Heskes T, Posthuma D. 2015.
  [MAGMA: generalized gene-set analysis of GWAS data](https://doi.org/10.1371/journal.pcbi.1004219).
- Mosquera JV, Auguste G, Wong D, et al. 2023.
  [Integrative single-cell meta-analysis reveals disease-relevant vascular cell states and markers in human atherosclerosis](https://doi.org/10.1016/j.celrep.2023.113380).
- Siewert A, Hoeland S, Mangold E, Ludwig KU. 2024.
  [Combining genetic and single-cell expression data reveals cell types and novel candidate genes for orofacial clefting](https://doi.org/10.1038/s41598-024-77724-9).
- Townsend HA, Rosenberger KJ, Vanderlinden LA, et al. 2024.
  [Evaluating methods for integrating single-cell data and genetics to understand inflammatory disease complexity](https://doi.org/10.3389/fimmu.2024.1454263).
- Mathys H, Boix CA, Akay LA, et al. 2024.
  [Single-cell multiregion dissection of Alzheimer's disease](https://doi.org/10.1038/s41586-024-07606-7).
- Lee D, Vicari JM, Porras C, et al. 2024 preprint.
  [Plasticity of Human Microglia and Brain Perivascular Macrophages in Aging and Alzheimer's Disease](https://doi.org/10.1101/2023.10.25.23297558).
- Zhang Y, et al. 2025.
  [Identifying Myeloid-Derived Suppressor Cells and Lipocalin-2 as Therapeutic Targets for Intervertebral Disc Degeneration](https://doi.org/10.1002/advs.202500505).
- Zhou M, Xue C. 2025.
  [Single-Cell Mapping of Genetic Risk Across Ten Respiratory Diseases](https://doi.org/10.3390/biology14121765).
- Mandelia M, Law PJ, Mills C, et al. 2025.
  [Deciphering genetic susceptibility to clear cell renal cell carcinoma](https://doi.org/10.1038/s42003-025-09297-w).
- Lai Q, Dannenfelser R, Roussarie JP, Yao V. 2025.
  [Disentangling associations between complex traits and cell types with seismic](https://doi.org/10.1038/s41467-025-63753-z).
- Li M, et al. 2025 preprint.
  [Benchmarking methods integrating GWAS and single-cell transcriptomic data for mapping trait-cell type associations](https://doi.org/10.1101/2025.05.24.25328275).
- Yao S, Harder A, Darki F, et al. 2025.
  [Connecting genomic results for psychiatric disorders to human brain cell types and regions reveals convergence with functional connectivity](https://doi.org/10.1038/s41467-024-55611-1).
- Li DY, Kundu S, Cheng P, et al. 2026.
  [Vascular smooth muscle cell state trajectories mediate molecular mechanisms of coronary disease risk](https://doi.org/10.1038/s41467-026-70530-z).
- [scDRS upstream repository and version notes](https://github.com/martinjzhang/scDRS).
- [scDRS CLI reference](https://martinjzhang.github.io/scDRS/reference_cli.html).
- [scDRS file formats](https://martinjzhang.github.io/scDRS/file_format.html).
- [CELLxGENE single-cell curation schema](https://chanzuckerberg.github.io/single-cell-curation/latest-schema.html).
- [CELLxGENE published-data download documentation](https://cellxgene.cziscience.com/docs/03__Download%20Published%20Data).
