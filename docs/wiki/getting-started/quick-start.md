# Quick Start

PostGWAS has two distinct entry points. Raw GWAS summary statistics first go
through the standalone `harmonisation` command. The downstream `pipeline`
starts from the resulting harmonised GWAS-VCF; it does not ingest the original
table directly.

## 1. Inspect the current interfaces

```console
postgwas harmonisation --help
postgwas pipeline --help
postgwas config --help
```

Do this with the exact checkout or installed version you will run.

## 2. Prepare a harmonisation sample sheet

Create one row per dataset. At minimum, map the input file, dataset identifier,
coordinates, effect and non-effect alleles, effect statistic, P value, allele
frequency, imputation quality, sample-size information, and trait type. Some
fields can be supplied as constants or external lookup files; the accepted
combinations are described in [Input Data](../reference/input-data.md).

Repository examples are available at:

```text
examples/configs/harmonisation/sample_sheet_quantitative.csv
examples/configs/harmonisation/sample_sheet_case_control.csv
examples/configs/harmonisation/run_config.yaml
```

Copy an example into your project and replace all study-specific values. Do not
run an example unchanged against unrelated data.

## 3. Validate, then harmonise

Use the exact syntax printed by `postgwas harmonisation --help`. The
harmonisation configuration and sample sheet are separate inputs. Resolve all
validation errors before starting a large run, then inspect the dataset-level
QC and rejection reports before using the GWAS-VCF downstream.

## 4. Configure downstream analyses

The pipeline configuration selects targets and supplies module settings. Check
the resolved configuration before execution:

```console
postgwas config validate --config /absolute/path/to/pipeline.yaml
postgwas config show --module filtering
```

The planner expands requested targets into their registered prerequisites. For
example, fine-mapping requires formatted inputs and LD clumping; the plan may
therefore contain stages that were not named directly.

## 5. Run from the harmonised GWAS-VCF

Follow the live `postgwas pipeline --help` interface to supply the harmonised
VCF, dataset identifier, output directory, resources, configuration, and target
modules. Before launching, verify genome build, ancestry, allele convention,
sample-size definition, and external-tool availability.

## 6. Review the run as a scientific record

Do not interpret only the final association table or plot. Retain the resolved
configuration, canonical log, commands, tool versions, stage summaries,
exclusion counts, and failure reports. See [Logging and Reproducibility](../core/logging-and-reproducibility.md)
and [Output Structure](../reference/output-structure.md).
