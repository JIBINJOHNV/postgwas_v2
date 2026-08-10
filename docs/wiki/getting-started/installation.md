# Installation

PostGWAS combines Python orchestration with external bioinformatics programs.
The repository's Dockerfile is the only installation definition that assembles
the complete tool stack. A Python or Conda installation is useful for
development and command inspection, but does not by itself provide every
external executable, R package, or reference dataset required by all modules.

## Recommended: build the container locally

From the repository root:

```console
docker build --platform linux/amd64 -t postgwas:local .
docker run --rm --platform linux/amd64 postgwas:local postgwas --help
```

The explicit platform matches the binaries installed by the current
Dockerfile. Building the image downloads licensed or third-party software;
review the Dockerfile and the terms for each tool before distributing an image.
In particular, do not assume the bundled MAGMA binary may be redistributed.

K-POPS and CALDERA are installed from immutable upstream commits into the
existing `postgwas` conda environment. K-POPS is installed as an environment
command, while CALDERA is installed below the environment's shared data
directory. PyTorch is added for K-POPS, and CALDERA uses the existing R,
`data.table`, and `dplyr` installation. No second K-POPS or CALDERA conda
environment is created.

Mount study data and resources rather than copying them into the image:

```console
docker run --rm --platform linux/amd64 \
  -v /absolute/path/to/project:/work \
  -v /absolute/path/to/resources:/resources:ro \
  postgwas:local \
  postgwas config validate --config /work/run.yaml
```

Use paths visible inside the container in YAML files and commands (`/work` and
`/resources` in this example), not host-only paths.

## Local development installation

Create the repository environment and install the checked-out package:

```console
conda env create -f environment.yml
conda activate postgwas
python -m pip install -e .
postgwas --help
```

This makes the Python CLI available. Before running a module, also install and
configure every external tool named in that module's guide. The `analysis` and
`enrichment` Python extras do not replace those system-level requirements.

For local K-POPS or CALDERA runs, prepare the commits recorded in the
Dockerfile with the repository's resource-preparation script. When supplied
the active environment's Python path, that script installs the pinned
`k-pops.py` command beside Python and CALDERA below the environment's shared
data directory. PostGWAS resolves both automatically, while its CALDERA R
adapter is resolved from the installed Python package.

## Verify the installation

Start with command discovery and configuration validation:

```console
postgwas --help
postgwas config --help
postgwas config validate --config /absolute/path/to/run.yaml
```

Validation checks configuration; it is not proof that all reference files or
external programs needed by a later analysis are scientifically compatible.
Run the relevant module's preflight checks before interpreting results.

## What is not installed automatically

Reference genomes, dbSNP files, allele-frequency panels, LD panels, gene
annotations, and other scientific resources are not interchangeable software
dependencies. They must be prepared separately and must match the configured
genome build and population. See [Reference Resources](../reference/reference-resources.md).
