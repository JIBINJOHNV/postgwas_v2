# PostGWAS

PostGWAS is a modular command-line software package for GWAS summary-statistics
harmonisation and downstream analysis. Harmonisation validates and standardises
variant coordinates, alleles and statistical fields, performs genome-build
conversion, and produces harmonised GRCh37 and GRCh38 GWAS-VCFs alongside
quality-control reports and records of transformations and variant exclusions.

Downstream analyses include variant filtering, summary-statistics imputation,
LD-block annotation, clumping, conditional and joint association analysis,
fine-mapping, gene and gene-set association testing, gene prioritisation, tissue
and single-cell analyses, pathway enrichment, SNP heritability estimation,
polygenic architecture modelling and visualisation.

Analyses can run as individual modules or supported downstream pipelines.
Pipeline mode requires harmonised summary statistics in PostGWAS GWAS-VCF format,
the necessary reference and resource files, and the selected analysis and run
settings. PostGWAS determines the required upstream steps, generates tool-specific
input files and executes the analyses in dependency order, automatically passing
outputs between stages. The resolved configuration, validation results, executed
commands and completion records are retained for reviewing and reproducing each
analysis.

[Install](#installation) · [Available analyses](#available-analyses) ·
[Start here](#choose-your-starting-point) ·
[Harmonisation](#harmonise-your-summary-statistics) ·
[Pipeline and direct mode](#run-downstream-analyses) ·
[Resources](#reference-data-and-configuration) ·
[Results](#results-logging-and-reproducibility) ·
[Full user guide](docs/wiki/home.md)

> Check genome build, ancestry, allele and gene identifiers, sample-size
> definitions, reference resources and QC before interpreting results. A
> successful command alone does not establish that an analysis is appropriate
> for your data.

## Installation

### Install PostGWAS and its analysis tools

The complete installer targets Apple-silicon macOS and glibc Linux x86-64.
Install and initialise [Miniforge with Mamba](https://github.com/conda-forge/miniforge)
first. macOS also requires Apple Command Line Tools and Rosetta 2.

```bash
git clone https://github.com/JIBINJOHNV/postgwas_v2.git
cd postgwas_v2
bash tools/setup/install_postgwas.sh --all-tools --name postgwas
conda activate postgwas
```

The installer creates a **new** environment. If `postgwas` already exists,
choose another name with `--name` and activate that name instead. You activate
one environment; tools with incompatible dependencies use managed internal
environments. Not every tool runs natively on Apple silicon.

Intel macOS supports the portable subset, not the complete installer. Native
Windows and Linux ARM64 are not supported by these installation paths.

### Verify the installation

The complete installer runs the repository's software verifier. These additional
checks confirm Python dependency consistency and public-command availability:

```bash
python -m pip check
postgwas --help
```

Software checks are not an analysis-readiness certificate. Reference panels,
expression data, annotation caches and service credentials must be prepared
separately for the analyses you select.

[Complete installation guide and platform table](docs/wiki/getting-started/installation.md) ·
[Resource setup](docs/wiki/getting-started/resource-setup.md) ·
[Troubleshooting](docs/wiki/help/troubleshooting.md)

## Available analyses

Every module below has a direct command. The **Pipeline target** column gives
the name accepted after `postgwas pipeline --modules`. Click a module for its
required inputs, resources, methods, commands, outputs and interpretation.

### Data preparation and quality control

| Module | Purpose | Pipeline target |
|---|---|---|
| [`harmonisation`](docs/wiki/harmonisation/overview.md) | Standardise raw summary statistics and produce harmonised GWAS-VCFs | Standalone only |
| [`sumstat_filter`](docs/wiki/modules/filtering.md) | Apply variant filters and record exclusions | `sumstat_filter` |
| [`qc`](docs/wiki/modules/qc-summary.md) | Assess VCF quality without creating a filtered VCF | `qc_summary` |
| [`formatter`](docs/wiki/modules/formatting.md) | Create MAGMA, GCTA, SuSiE, FINEMAP, PRED-LD, LDSC and MiXeR inputs | `formatter` |
| [`imputation`](docs/wiki/modules/imputation.md) | Impute summary statistics using PRED-LD/TOP_LD and re-harmonise the results | `imputation` |

### Locus analysis and fine-mapping

| Module | Available analyses | Pipeline target |
|---|---|---|
| [`annot_ldblock`](docs/wiki/modules/ld-annotation.md) | LD-block annotation | `annot_ldblock` |
| [`ld_clump`](docs/wiki/modules/ld-clumping.md) | Region-based, standard and COJO-selection clumping | `ld_clump` |
| [`gcta_cojo`](docs/modules/gcta_cojo/README.md) | Stepwise selection, top-SNP, joint and conditional analysis | `gcta_cojo` |
| [`finemap`](docs/wiki/modules/fine-mapping.md) | SuSiE-RSS and FINEMAP | `finemap` |

### Gene, gene-set and cell-type analysis

| Module | Available analyses | Pipeline target |
|---|---|---|
| [`magma`](docs/wiki/modules/magma.md) | Positional MAGMA, eMAGMA, H-MAGMA, nMAGMA and chromMAGMA mappings; optional gene-set tests | `magma` |
| [`gcta_gene`](docs/wiki/modules/gcta-gene.md) | fastBAT gene, segment and set tests; mBAT-combo | `gcta_gene` |
| [`magmacovar`](docs/wiki/modules/magmacovar.md) | MAGMA gene-property and conditional analyses | `magmacovar` |
| [`single_cell`](docs/wiki/modules/single-cell.md) | MAGMA cell-type analysis, scDRS and LDSC cell-type analysis | `single_cell` |

### Gene prioritisation

| Module | Approach | Pipeline target |
|---|---|---|
| [`pops`](docs/wiki/modules/pops.md) | Feature-based gene prioritisation | `pops` |
| [`kpops`](docs/modules/kpops.md) | Kernel-based gene prioritisation | `kpops` |
| [`caldera`](docs/modules/caldera.md) | Integrate PoPS predictions and credible sets | `caldera` |
| [`flames`](docs/wiki/modules/flames.md) | Integrate fine-mapping, MAGMA, gene-property and PoPS evidence | `flames` |

### Heritability, polygenicity and interpretation

| Module | Available analyses | Pipeline target |
|---|---|---|
| [`heritability`](docs/wiki/modules/ldsc.md) | Single-trait LDSC heritability, with optional liability conversion | `heritability` |
| [`mixer`](docs/wiki/modules/mixer.md) | Univariate MiXeR, GSA-MiXeR, or both | `mixer` |
| [`pathway_enrichment`](docs/wiki/modules/pathway-enrichment.md) | Gene-list enrichment and interaction services | Standalone only |
| [`manhattan`](docs/wiki/modules/manhattan.md) | Manhattan plots | `manhattan` |

Harmonisation and pathway enrichment are **standalone-only**. `qc` is the direct
command and `qc_summary` is its pipeline target. Supporting interfaces are
`postgwas config`, `postgwas resources`, `postgwas pipeline`, and
`postgwas --validate`; see the [command reference](docs/wiki/reference/command-reference.md).

## Choose your starting point

### I have raw GWAS summary statistics

Prepare a sample sheet describing your files, columns and study metadata, then
run harmonisation. Review the resulting QC and rejection evidence before
selecting downstream analyses.

[Sample-sheet guide](docs/wiki/harmonisation/sample-sheet.md) ·
[Connected quick-start tutorial](docs/wiki/getting-started/quick-start.md)

### I already have a PostGWAS-harmonised VCF

Use pipeline mode to request downstream analyses, or run a direct module that
accepts VCF input. Pipeline entry requires an indexed, single-sample
PostGWAS-harmonised VCF with the required provenance; an arbitrary GWAS-VCF is
not sufficient. Do not add provenance headers manually to bypass this boundary.

[Input contracts](docs/wiki/core/input-output-contracts.md) ·
[Pipeline input validation](docs/wiki/core/pipeline-input-validation.md)

### I already have tool-specific analysis files

Use the appropriate direct command—for example, MAGMA results for PoPS, or
formatted statistics and a locus table for fine-mapping. Direct inputs are not
universally VCFs. Check the selected module's build, identifier and resource
requirements even when another program produced those inputs.

[Running modules independently](docs/wiki/core/running-modules-independently.md)

## Harmonise your summary statistics

### Prepare the sample sheet and resources

A version-2 CSV/TSV sample sheet maps each study's columns and metadata to
PostGWAS concepts. It is separate from the run-configuration YAML, which
controls processing policies and resources. Complete the required EAF, INFO,
effect and sample-size source choices before running; generated drafts still
require your review.

[Sample-sheet requirements](docs/wiki/harmonisation/sample-sheet.md) ·
[Configuration guide](docs/modules/harmonisation/configuration.md) ·
[Complete source, step and policy reference](docs/wiki/harmonisation/policies.md)

### Run harmonisation

The connected example uses dataset `STUDY`. Replace `studies.csv` with your
completed sample sheet containing that identifier, and prepare the referenced
resources first. Omit `--dataset-id` to process all sample-sheet rows.

```bash
postgwas harmonisation \
  --sample-sheet studies.csv \
  --dataset-id STUDY \
  --resource-directory reference/harmonisation \
  --output-directory results/harmonisation \
  --validate
```

Harmonisation prepares each dataset once, processes its chromosomes, then
merges and assesses the outputs. Study-wide decisions—including build, strand
consensus, effect type, OR SE scale, P-value scale and EAF/MAF interpretation—are
recorded as `DECIDE` entries in the log. Dataset rows currently run sequentially;
chromosome work uses bounded parallelism.

### Review the results before continuing

Successful runs produce GRCh37 and GRCh38 VCFs with indexes, QC reports and
rejection evidence. For this example, the GRCh37 file is
`results/harmonisation/STUDY/harmonisation/STUDY_GRCh37_merged.vcf.gz`.

Harmonisation can transform or reject variants, including unresolved
palindromes and unusable effect statistics. Inspect `rejected/`, the
reason-by-chromosome matrix and liftover loss counts. The post-merge QC
assessment reports without removing additional variants; `sumstat_filter`
materialises later QC-based exclusions.

`--validate` adds an original-input versus same-build output concordance check.
Review that result and the dataset/run completion status before continuing.

[All processing steps](docs/wiki/harmonisation/processing-order.md) ·
[Outputs and QC](docs/wiki/harmonisation/outputs-and-qc.md) ·
[Concordance validation](docs/wiki/reference/validation.md)

## Run downstream analyses

### Pipeline mode: start from a harmonised VCF

This positional MAGMA example assumes a GRCh37 study, a compatible EUR PLINK
reference whose BIM uses rsIDs, and NCBI37.3 gene locations with Entrez
identifiers. These are example resource choices, not suitable defaults for
every study. The prefix must identify matching `.bed`, `.bim` and `.fam` files.

```bash
postgwas pipeline \
  --modules magma \
  --vcf results/harmonisation/STUDY/harmonisation/STUDY_GRCh37_merged.vcf.gz \
  --magma-ld-reference reference/GRCh37_EUR_reference \
  --gene-location-file reference/NCBI37.3.gene.loc \
  --resolve-variants-to-reference \
  --dataset-id STUDY \
  --output-directory results/magma_pipeline
```

PostGWAS creates the MAGMA input tables using the validated BIM identifier
convention and passes them to MAGMA. The explicit variant-resolution option
retains exact identifier matches to the BIM; it is not allele harmonisation or
liftover. Inspect its overlap and exclusion audit. The pipeline does not
download the reference genotypes or gene locations.

### Direct mode: start from the module's required files

For the same analysis, first create MAGMA input tables if you do not already
have them:

```bash
postgwas formatter \
  --vcf results/harmonisation/STUDY/harmonisation/STUDY_GRCh37_merged.vcf.gz \
  --format magma \
  --variant-id-type rsid \
  --dataset-id STUDY \
  --output-directory results/formatted
```

Then run MAGMA directly with those paired SNP-location and p-value/sample-size
tables and the same references:

```bash
postgwas magma \
  --snp-location-file results/formatted/STUDY_magma_snp_loc.tsv \
  --p-value-file results/formatted/STUDY_magma_p_values.tsv \
  --magma-ld-reference reference/GRCh37_EUR_reference \
  --gene-location-file reference/NCBI37.3.gene.loc \
  --resolve-variants-to-reference \
  --dataset-id STUDY \
  --output-directory results/magma_direct
```

Direct mode runs the selected module's documented workflow and validations;
some modules still prepare internal intermediates. Pipeline mode additionally
plans registered upstream modules and passes their artifacts between stages.

[MAGMA methods and interpretation](docs/wiki/modules/magma.md) ·
[Direct-mode guide](docs/wiki/core/running-modules-independently.md)

### Inspect the selected workflow

```bash
postgwas magma --help
postgwas pipeline --modules magma --help
```

Dependencies change with the selected method: region clumping needs LD-block
annotation, standard clumping uses prepared pairwise LD, and COJO selection
needs formatted statistics and PLINK genotypes. LDSC cell typing does not use
the MAGMA branch. Fine-mapping and gene-prioritisation targets can combine
several branches; selecting one does not mean every module runs.

[Pipeline dependencies and execution order](docs/wiki/core/pipeline-workflow.md) ·
[Single-cell methods](docs/wiki/modules/single-cell.md)

## Reference data and configuration

### Prepare compatible reference resources

PLINK genotypes, indexed pairwise-LD tables, LD-block intervals, LDSC scores,
expression matrices and K-POPS kernels are distinct resource types. Match the
genome build, ancestry, allele conventions and gene-identifier system required
by each selected method. A pipeline creates intermediate analysis files, not
all of these external resources.

A complete universal reference bundle is not currently published. Follow the
module-specific preparation instructions. One supported downloader installs
the pinned MAGMA functional-mapping bundle and writes its configuration
fragment; it does not replace other required MAGMA or downstream references:

```bash
postgwas resources prepare magma --output-directory reference/magma_mappings
```

[Resource setup](docs/wiki/getting-started/resource-setup.md) ·
[Reference requirements](docs/wiki/reference/reference-resources.md)

### Use CLI options or a run configuration

Commands supporting `--run-config` resolve packaged defaults, user YAML and
explicit CLI overrides in that order. Export a template, fill in your resource
paths, and validate it before analysis:

```bash
postgwas config export --pipeline magma --output magma_pipeline.yaml
postgwas config validate --config magma_pipeline.yaml
```

Configuration validation checks the schema, not the suitability of all input
data. Direct `manhattan` and `pathway_enrichment` do not expose `--run-config`;
consult their CLI options. Configuration module names can differ from public
commands, such as `formatting` versus `formatter`.

[Configuration guide](docs/wiki/core/configuration.md) ·
[Packaged settings](docs/wiki/reference/configuration-defaults.md)

## Results, logging and reproducibility

### Find and interpret your results

Harmonisation groups files by dataset. Pipelines create numbered step
directories under the requested output directory; the numbers follow the
actual plan rather than fixed module identifiers. Native outputs, prepared
inputs, logs and normalised results follow each module's documented layout.

Outputs include VCFs, association tables, credible sets, gene scores, images
and native-tool files. HTML reports are available for some modules, not all.
Association p-values, PoPS/K-POPS ranking scores and CALDERA probabilities have
different interpretations; chromMAGMA rankings are not calibrated gene
association p-values.

[Output structure and archiving](docs/wiki/reference/output-structure.md)

### Monitor, resume or restart an analysis

Terminal display is on by default. `--hide-screen` hides the terminal copy but
retains the saved screen transcript. Keep the resolved configuration, commands,
software/resource versions, QC evidence and native-tool logs with your results.

`--resume` is the default: existing results are reused only after checkpoint
validation. Incomplete work restarts at a safe boundary, and changed inputs or
parameters trigger the configured warning-and-restart policy. `--overwrite`
requests a restart from the beginning. Replacement remains restricted to
unchanged, recorded PostGWAS-owned outputs; modified or unrecognised files are
preserved and cause an actionable stop.

[Logging and reproducibility](docs/wiki/core/logging-and-reproducibility.md) ·
[Resume and overwrite policy](docs/wiki/core/configuration.md#resume-and-overwrite)

## Documentation and support

### Troubleshooting and current limitations

Start with the first validation or tool error and inspect the relevant QC
report before changing thresholds. Important current boundaries include:

- Harmonisation and pathway enrichment are standalone-only. The pipeline
  accepts one PostGWAS-harmonised study VCF per invocation.
- Differently anchored or padded indels are not FASTA-normalised before the
  strand-reference join and may be rejected as unmatched. See the
  [harmonisation processing guide](docs/wiki/harmonisation/processing-order.md).
- CALDERA pipeline integration currently supports GRCh37 only. Direct-mode
  GRCh38 has additional identifier requirements described in its guide.
- Heritability and MiXeR currently expose single-trait analyses, not genetic
  correlation or bivariate MiXeR. The Manhattan module does not generate QQ plots.
- Some clumping/fine-mapping combinations with `--apply-filter` currently fail
  during argument parsing because of a duplicate `--remove-mhc` option.
  This is not a general prohibition on filtering all pipelines. See the
  [specific limitation and separate-filter workaround](docs/wiki/core/pipeline-workflow.md).
- Enrichment providers require network access and sometimes credentials;
  inspect each provider's outcome rather than treating command success as
  proof that all providers succeeded.

[Troubleshooting](docs/wiki/help/troubleshooting.md) ·
[Error messages](docs/wiki/help/error-messages.md) ·
[Frequently asked questions](docs/wiki/help/faq.md)

### Full documentation and method references

The [PostGWAS user guide](docs/wiki/home.md) provides the complete documentation
index, including every module and the detailed harmonisation policy reference.
The pages are version-controlled Markdown and can be browsed when GitHub Wiki
access is not enabled. Each module guide explains inputs, direct/pipeline
commands, outputs, assumptions and limitations.

Cite the methods, software and reference datasets used in your analysis.

[Method references](docs/wiki/reference/scientific-references.md) ·
[Analysis assumptions and limitations](docs/wiki/core/scientific-considerations.md)

## Development and licensing

### Documentation and development checks

```bash
python -m pytest -q
python tools/docs/build_wiki.py --check
python tools/docs/validate_wiki_cli.py
python tools/docs/update_harmonisation_policies.py --check
```

Edit canonical documentation sources, not generated Wiki copies. Navigation
comes from `docs/wiki.yml`; the full index lives in the user-guide home page.
Repository changes require method validation, regression tests and protection
of the bundled harmonisation adapters.

[Documentation maintenance](docs/wiki/README.md) ·
[Module-page template](docs/templates/module-page.md)

### Licensing status

PostGWAS does not currently publish a repository licence. Until the copyright
holder selects and adds one, do not assume permission to redistribute or modify
the project. Third-party tools and reference datasets carry their own terms,
which apply independently.
