# PostGWAS

PostGWAS is a command-line toolkit for converting raw GWAS summary statistics
into harmonised GWAS-VCF and running reproducible post-GWAS analyses. Its
modules cover quality control, imputation, LD-based analyses, fine-mapping,
heritability, gene and pathway analysis, gene prioritisation, and single-cell
integration.

> PostGWAS is research software. Review the resolved configuration, scientific
> assumptions, QC reports, rejection evidence, and external-tool output before
> interpreting results.

## Documentation

The complete user guide is stored as version-controlled Markdown in this
repository, so it remains available when GitHub Wiki access is not enabled.
Start with the [PostGWAS User Guide](docs/wiki/home.md), or open a topic below.

### Getting started

- [Home](docs/wiki/home.md)
- [Installation](docs/wiki/getting-started/installation.md)
- [Quick Start](docs/wiki/getting-started/quick-start.md)
- [Resource Setup](docs/wiki/getting-started/resource-setup.md)

### Core concepts

- [Configuration](docs/wiki/core/configuration.md)
- [Pipeline Workflow](docs/wiki/core/pipeline-workflow.md)
- [Input and Output Contracts](docs/wiki/core/input-output-contracts.md)
- [Logging and Reproducibility](docs/wiki/core/logging-and-reproducibility.md)
- [Scientific Considerations](docs/wiki/core/scientific-considerations.md)
- [Running Modules Independently](docs/wiki/core/running-modules-independently.md)

### Harmonisation

- [Harmonisation Overview](docs/wiki/harmonisation/overview.md)
- [Harmonisation Sample Sheet](docs/wiki/harmonisation/sample-sheet.md)
- [Harmonisation Configuration](docs/modules/harmonisation/configuration.md)
- [Harmonisation Processing Order](docs/wiki/harmonisation/processing-order.md)
- [Harmonisation Outputs and QC](docs/wiki/harmonisation/outputs-and-qc.md)

### Analysis modules

- [Filtering](docs/wiki/modules/filtering.md)
- [Formatting](docs/wiki/modules/formatting.md)
- [QC Summary](docs/wiki/modules/qc-summary.md)
- [Manhattan Plots](docs/wiki/modules/manhattan.md)
- [Imputation](docs/wiki/modules/imputation.md)
- [LD Annotation](docs/wiki/modules/ld-annotation.md)
- [LD Clumping](docs/wiki/modules/ld-clumping.md)
- [LDSC Heritability](docs/wiki/modules/ldsc.md)
- [Fine Mapping](docs/wiki/modules/fine-mapping.md)
- [MAGMA](docs/wiki/modules/magma.md)
- [GCTA Gene Analysis](docs/wiki/modules/gcta-gene.md)
- [MAGMAcovar](docs/wiki/modules/magmacovar.md)
- [Single-Cell Integration](docs/wiki/modules/single-cell.md)
- [PoPS](docs/wiki/modules/pops.md)
- [K-POPS](docs/modules/kpops.md)
- [CALDERA](docs/modules/caldera.md)
- [FLAMES](docs/wiki/modules/flames.md)
- [MiXeR](docs/wiki/modules/mixer.md)
- [Pathway Enrichment](docs/wiki/modules/pathway-enrichment.md)

### Reference

- [Input Data](docs/wiki/reference/input-data.md)
- [Reference Resources](docs/wiki/reference/reference-resources.md)
- [Configuration Defaults](docs/wiki/reference/configuration-defaults.md)
- [Output Structure](docs/wiki/reference/output-structure.md)
- [Command Reference](docs/wiki/reference/command-reference.md)
- [Validation Reference](docs/wiki/reference/validation.md)
- [Scientific References](docs/wiki/reference/scientific-references.md)
- [scDRS Evidence Review](docs/wiki/reference/scdrs-evidence-review.md)

### Help

- [Troubleshooting](docs/wiki/help/troubleshooting.md)
- [Frequently Asked Questions](docs/wiki/help/faq.md)
- [Error Messages](docs/wiki/help/error-messages.md)

## Installation

The repository container is the complete installation definition. Build it
locally from the repository root:

```console
docker build --platform linux/amd64 -t postgwas:local .
docker run --rm --platform linux/amd64 postgwas:local postgwas --help
```

For Python development and command inspection:

```console
conda env create -f environment.yml
conda activate postgwas
python -m pip install -e .
postgwas --help
```

A local Python installation does not install every external bioinformatics
program or scientific reference resource. See the
[installation guide](docs/wiki/getting-started/installation.md) and
[resource setup guide](docs/wiki/getting-started/resource-setup.md).

## Quick start

Raw summary statistics first pass through the standalone harmonisation command.
Start from one of the maintained version-2 sample-sheet templates:

```text
examples/configs/harmonisation/sample_sheet_quantitative.csv
examples/configs/harmonisation/sample_sheet_case_control.csv
examples/configs/harmonisation/run_config.yaml
```

Then inspect and run the live interface:

```console
postgwas harmonisation --help

postgwas harmonisation \
  --sample-sheet studies.csv \
  --run-config harmonisation.yaml \
  --resource-directory /absolute/path/to/resources \
  --output-directory /absolute/path/to/results
```

Downstream workflows begin from a successfully harmonised GWAS-VCF:

```console
postgwas pipeline --help
postgwas config --help
```

PostGWAS resolves requested pipeline targets and their registered dependencies.
Do not copy old command examples or infer option names; use the help output from
the installed version and the [Quick Start](docs/wiki/getting-started/quick-start.md).

## Configuration and reproducibility

Packaged YAML provides canonical defaults. User YAML overrides packaged YAML,
and explicit CLI values override user YAML. Export or inspect current settings
instead of copying defaults from older runs:

```console
postgwas config export \
  --module harmonisation \
  --style full \
  --output harmonisation.yaml

postgwas config validate --config harmonisation.yaml
```

Keep the resolved configuration, run manifest, canonical logs, tool versions,
commands, QC summaries, rejection reports, and reference-resource provenance
with every analysis.

## Development checks

```console
python -m pytest -q
python tools/docs/build_wiki.py --check
python tools/docs/validate_wiki_cli.py
```

Repository development must follow [AGENTS.md](AGENTS.md), including scientific
validation, focused regression tests, cumulative-diff review, and protection of
the bundled harmonisation adapters.
