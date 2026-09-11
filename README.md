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

## Workflow at a glance

```text
Raw summary statistics + reviewed sample sheet
  → standalone harmonisation → GRCh37 and GRCh38 GWAS-VCFs + QC
                                 ↓ choose one compatible build
                         Selected downstream pipeline
                           ├─ Variant and locus analyses
                           ├─ Gene, gene-set and cell-type analyses
                           └─ Heritability, polygenicity and plots
```

You supply the study data and the selected methods' reference/resource files.
The pipeline prepares tool-specific inputs and connects the required analysis
steps. Select one or more branches; not every module needs to run.
Compatible existing tool-specific inputs can instead enter a direct module;
gene-list pathway enrichment is also a separate direct workflow.

> Check genome build, ancestry, allele and gene identifiers, sample-size
> definitions, reference resources and QC before interpreting results. A
> successful command alone does not establish that an analysis is appropriate
> for your data.

## Installation

Choose a local Mamba installation or build the Linux Docker image. Both routes
require the reference data for your selected analyses separately.

### Option 1: install locally with Mamba

The full installer installs PostGWAS and its analysis tools on Apple-silicon
Macs and glibc-based x86-64 Linux systems. Before running it, install and
initialise [Miniforge with Mamba](https://github.com/conda-forge/miniforge).
Apple-silicon Macs also require Apple Command Line Tools and Rosetta 2.

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

### Verify the local installation

The complete installer runs the repository's software verifier. These additional
checks confirm Python dependency consistency and public-command availability:

```bash
python -m pip check
postgwas --help
```

Software checks are not an analysis-readiness certificate. Reference panels,
expression data, annotation caches and service credentials must be prepared
separately for the analyses you select.

### Option 2: build and run with Docker

Install and start [Docker Desktop on macOS](https://docs.docker.com/desktop/setup/install/mac-install/)
or [Docker Engine on Linux](https://docs.docker.com/engine/install/). From a
clean checkout of this repository, build the supplied [Dockerfile](Dockerfile):

```bash
docker build --platform linux/amd64 --load --tag postgwas:local .
docker run --rm --platform linux/amd64 postgwas:local postgwas --help
```

No host Mamba environment is needed for this route. The image targets Linux
x86-64; Apple-silicon Macs require amd64 emulation, which can be slower. It is
not a native ARM64 image. See [Docker's platform guidance](https://docs.docker.com/build/building/multi-platform/).
The Dockerfile runs software checks during the build; these are not end-to-end
analysis tests, and current CI does not build or test the Docker image.

Mount input/reference folders and a writable output folder when running an
analysis. Paths in commands, sample sheets and YAML must refer to locations
inside the container. The guide explains persistent results, permissions,
build-context privacy and third-party image redistribution restrictions:

[Docker setup and running analyses](docs/wiki/getting-started/installation.md#docker-installation)

[Complete installation guide and platform table](docs/wiki/getting-started/installation.md) ·
[Resource setup](docs/wiki/getting-started/resource-setup.md) ·
[Troubleshooting](docs/wiki/help/troubleshooting.md)

### Before your first analysis

- **Software:** complete one installation route and its software checks.
- **Study inputs:** keep the original statistics and study documentation; review
  a sample sheet for raw inputs, or check the existing VCF's provenance and index.
- **References:** prepare the build-, population- and method-matched resources
  listed in the selected module guides. They are separate from the software.
- **Storage and compute:** use a writable output location with room for results,
  intermediates and logs, and fit thread/memory settings to the host or container.

The [resource checklist](docs/wiki/getting-started/resource-setup.md#preflight-checklist)
and [connected tutorial](docs/wiki/getting-started/quick-start.md) take you from
installation to a first analysis using your own data.

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

### Quick start: from raw statistics to results

The connected example uses one dataset, `STUDY`, throughout:

1. [Prepare the inputs and sample sheet](#generate-and-review-the-sample-sheet),
   including the study's real sample sizes and compatible references.
2. [Run harmonisation](#run-harmonisation) and review both build outputs,
   rejected variants and concordance evidence.
3. [Assess the resulting VCF with QC](docs/wiki/getting-started/quick-start.md#4-review-qc-without-changing-the-vcf).
   QC creates an assessment, not a filtered VCF.
4. [Run positional MAGMA in pipeline mode](#pipeline-mode-start-from-a-harmonised-vcf).
   The pipeline prepares the MAGMA inputs; review the gene results and exclusions.

Follow the [full walkthrough](docs/wiki/getting-started/quick-start.md) for the
connected commands and output paths. Its paths are placeholders for your study
and references, not a bundled demonstration dataset.

## Harmonise your summary statistics

### Generate and review the sample sheet

A version-2 CSV/TSV sample sheet maps each study's columns and metadata to
PostGWAS concepts. It is separate from the run-configuration YAML, which
controls processing policies and resources. Generate a draft from the headers
of files in your raw-data directory:

```bash
python -m postgwas.modules.harmonisation.sample_sheet_generator \
  --input-directory raw_data \
  --output studies.csv
```

Review every generated row against the original study documentation. The
generator can leave fields unresolved or omit files it cannot map; inspect its
warnings and `studies.csv.rejected_files.tsv` report. Successful draft generation
does not mean the studies are ready to run. Alternatively, copy a maintained
[quantitative template](examples/configs/harmonisation/sample_sheet_quantitative.csv)
or [case-control template](examples/configs/harmonisation/sample_sheet_case_control.csv)
and replace all placeholder values.

### Complete the study information

| Information | What to provide or check |
|---|---|
| Identity and file | `config_version` set to `2`, a unique `dataset_id`, and `input_file`; relative file paths resolve from the sample-sheet directory. |
| Coordinates and alleles | Chromosome/position columns (or a combined column), effect allele and other allele. Check that the effect and frequency refer to the declared effect allele. |
| Effect and P value | An effect estimate or Z score, a P-value column and its representation (`raw`, `neglog10` or `auto`). Map supplied SE when available; recovery of missing statistics is conditional. |
| Allele frequency | Exactly one internal `effect_allele_frequency_column` or external file/column pair. The default comparison AF panel does not replace this study-frequency source. |
| Imputation quality | Internal INFO takes priority over an external INFO file/column pair. With neither source, explicitly choose `--fixed-info VALUE`; an external proxy or fixed value is not study-measured quality. |
| Sample size | For quantitative traits, `control_count_column` or `control_count` represents total N. For case-control traits, provide real control and case counts through their corresponding column or fixed-count fields. |

Do not invent missing counts, frequencies or quality measurements. External
per-chromosome tables need explicit `{chromosome}` path templates, not just a
filename prefix. Prepare the harmonisation reference tree before running.

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

### Pipeline targets and execution order

Select the final analyses you want; PostGWAS adds their registered prerequisites.
All public pipeline targets are listed below. Each row shows a single-target
plan without optional workflow stages; arrows show scheduled order, not a claim
that each stage consumes only the immediately preceding stage's output.

For the `finemap`, `caldera` and `flames` rows, explicitly select
`--clumping-methods standard` and either `--finemap-method susie` or
`--finemap-method finemap`. The engines use the same stage order but different
tool-specific inputs and resources.

| Selected target and method | Execution order |
|---|---|
| `sumstat_filter` | `sumstat_filter` |
| `qc_summary` | `qc_summary` |
| `formatter` | `formatter` |
| `imputation` | `formatter → imputation` |
| `annot_ldblock` | `annot_ldblock` |
| `ld_clump`, `standard` only | `ld_clump` |
| `ld_clump`, `region` | `annot_ldblock → ld_clump` |
| `ld_clump`, `cojo-slct` | `formatter → ld_clump` |
| `finemap`, standard clumping | `ld_clump → formatter → finemap` |
| `magma` | `formatter → magma` |
| `gcta_gene` | `formatter → gcta_gene` |
| `gcta_cojo` | `formatter → gcta_cojo` |
| `magmacovar` | `formatter → magma → magmacovar` |
| `single_cell`, `magma_celltype` or `scdrs` | `formatter → magma → single_cell` |
| `single_cell`, `ldsc_celltype` only | `formatter → single_cell` |
| `pops` | `formatter → magma → pops` |
| `kpops` | `formatter → magma → kpops` |
| `caldera`, standard clumping | `formatter → ld_clump → magma → pops → finemap → caldera` |
| `flames`, standard clumping | `ld_clump → formatter → magma → magmacovar → pops → finemap → flames` |
| `heritability` | `formatter → heritability` |
| `mixer` | `formatter → mixer` |
| `manhattan` | `manhattan` |

FLAMES combines fine-mapping, MAGMA, MAGMAcovar and PoPS evidence. CALDERA
combines PoPS with fine-mapped credible sets; it does not run MAGMAcovar or
K-POPS, and its pipeline currently requires GRCh37. Harmonisation and pathway
enrichment remain standalone-only, and
`qc_summary` is the pipeline target for the direct `qc` command.

The packaged clumping selection is `region` plus `standard`, which also
schedules `annot_ldblock`. The standard-only rows above deliberately omit that
stage. Fine-mapping needs the standard-clumping locus artifact; adding region
or COJO analyses can change prerequisites and order, but neither replaces that
artifact. Always inspect the plan for your actual method combination.

For multiple targets, use the combined plan rather than concatenating table
rows. Shared prerequisites normally run once. Formatting before and after imputation
uses different data and can therefore appear twice. Inspect contextual help for
the methods you intend to use:

```bash
postgwas magma --help
postgwas pipeline --modules magma --help
postgwas pipeline --modules flames --clumping-methods standard --finemap-method susie --help
```

Dependencies change with the selected method: region clumping needs LD-block
annotation, standard clumping uses prepared pairwise LD, and COJO selection
needs formatted statistics and PLINK genotypes. LDSC cell typing does not use
the MAGMA branch. Fine-mapping and gene-prioritisation targets can combine
several branches; selecting one does not mean every module runs.

[Pipeline dependencies and execution order](docs/wiki/core/pipeline-workflow.md) ·
[Single-cell methods](docs/wiki/modules/single-cell.md)

### Optional pipeline stages

| Switch | Stage added |
|---|---|
| `--apply-filter` | Variant filtering; when combined with imputation, also a separate post-imputation filter. |
| `--apply-imputation` | Input formatting, PRED-LD imputation and internal re-harmonisation; requires PRED-LD and harmonisation resources. |
| `--apply-manhattan` | Manhattan plots for the pipeline's resulting study VCF. |
| `--heritability` | LDSC heritability estimation; requires the merge-allele file, LD scores and weights. |

These switches add work and may add resource requirements. They are not needed
for every analysis. MAGMA and some clumping/fine-mapping combinations currently
conflict with `--apply-filter`; read the
[supported combinations and separate-filter workaround](docs/wiki/core/pipeline-workflow.md#optional-workflow-stages)
before combining switches.

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

### Check reference compatibility

- Match the study and reference genome builds, chromosome labels and coordinate
  conventions; a filename alone is not evidence of compatibility.
- Match the LD/frequency reference population to the analysis and study design.
- Check variant IDs and allele conventions across summary statistics and LD
  resources, and gene IDs across gene locations, gene sets and expression data.
- Keep required indexes and companion files together; record reference sources,
  releases and checksums. Review overlap and exclusion counts after each join.

Use the [module-by-module reference checklist](docs/wiki/reference/reference-resources.md)
and [input contracts](docs/wiki/core/input-output-contracts.md) for the exact
requirements. Validation does not make an unsuitable reference appropriate.

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

### Check that the run succeeded

1. Confirm successful completion in the canonical log and completion records,
   including every requested pipeline stage. A file's existence is not enough.
2. Check the module's primary results and required companions, then open its
   HTML report where available. Compare input, retained and excluded counts.
3. Review warnings, native-tool diagnostics and join/allele checks. For
   harmonisation, inspect rejected records, the reason matrix and liftover losses.
4. Retain that evidence with the resolved configuration, commands and
   software/reference versions. Review whether the method and references fit
   the study before interpreting associations or rankings.

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
- MAGMA and some clumping/fine-mapping combinations with `--apply-filter`
  currently fail during argument parsing because of duplicate MHC options.
  This is not a general prohibition on filtering all pipelines. See the
  [specific limitation and separate-filter workaround](docs/wiki/core/pipeline-workflow.md).
- Enrichment providers require network access and sometimes credentials;
  inspect each provider's outcome rather than treating command success as
  proof that all providers succeeded.

[Troubleshooting](docs/wiki/help/troubleshooting.md) ·
[Error messages](docs/wiki/help/error-messages.md) ·
[Frequently asked questions](docs/wiki/help/faq.md)

### Moving from the previous PostGWAS version

Do not reuse an older installation command, sample sheet, YAML or result path
without checking the current interface. The
[migration guide](docs/wiki/getting-started/migrating-from-v1.md) maps the main
changes and explains how to rebuild a reviewed configuration. Keep the old run
separate; old outputs are not automatically valid version-2 checkpoints.

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
