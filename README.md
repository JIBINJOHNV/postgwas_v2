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

- [GitHub Wiki](https://github.com/JIBINJOHNV/postgwas_v2/wiki) — complete user
  guide, workflows, module pages, inputs, outputs, and troubleshooting.
- [Canonical documentation sources](docs/wiki/home.md) — the version-controlled
  source used to generate the Wiki.
- [Harmonisation guide](docs/wiki/harmonisation/overview.md) — preparing raw
  summary statistics and creating GWAS-VCF.
- [Live command reference](docs/wiki/reference/command-reference.md) — use with
  the `--help` output from the exact checkout being run.

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
