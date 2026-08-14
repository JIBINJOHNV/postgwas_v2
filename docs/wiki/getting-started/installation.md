# Installation

PostGWAS has three separate installation layers:

1. the PostGWAS Python package and its Python/R libraries;
2. third-party command-line programs used by selected modules; and
3. genome-build-, population-, and method-specific scientific resources.

Creating the Conda environment completes the first layer and installs the
portable command-line tools available for the host platform. It does **not**
make every scientific analysis ready. A module is ready only after every item
in its row of the [module requirements table](#module-requirements) is present
and compatible.

## Supported installation paths

| Host | Local environment | Complete external-tool stack |
|---|---|---|
| Linux x86-64 | Supported by `environment.yml`; Linux selectors add PLINK 2, GCTA, and BGENIX | Install the remaining licensed/upstream tools or build the container |
| macOS Apple Silicon | Tested for the Python/R stack, bcftools, tabix, bedtools, pigz, PLINK 1.9, scDRS, and the bcftools liftover plugin | Several upstream tools are Linux/x86-64 only; use the Linux container for workflows that require them |
| macOS Intel | Dependency solve verified; native runtime is not yet tested | Prefer the Linux container for Linux-only workflows |
| Windows | Native installation is not supported | Use Linux through WSL2 or a Linux container host |

The environment uses [CEP-24 platform selectors](https://conda.org/learn/ceps/cep-0024/).
Use a current Mamba or Conda release that understands dictionary selectors;
older clients may reject entries such as `sel(linux)`.

## Prerequisites

- Git and a complete checkout of this repository.
- Miniforge, Mambaforge, or another current Conda distribution with access to
  `conda-forge` and `bioconda`.
- At least 5 GB free for the local software environment. A tested Apple Silicon
  installation occupied about 2.5 GB before scientific reference data.
- Substantially more storage for reference genomes, LD panels, and study data.
- Docker with an active daemon when using the Linux container path.

Do not use the system Python. The all-module Python/R environment is pinned to
Python 3.10 because the current scanpy/scDRS/Numba stack is not compatible with
Python 3.12 or later.

## Install with one command

The recommended installer creates the Conda environment, installs PostGWAS,
builds the pinned bcftools liftover plugin inside the same environment, and
checks the resulting command-line installation. Run it from the repository
root:

```console
bash tools/setup/install_postgwas.sh
```

The script uses Mamba when available and otherwise uses Conda. It supports
Apple Silicon macOS, Intel macOS, and glibc-based Linux x86-64 distributions,
including Ubuntu, Debian, Rocky Linux, and Fedora. Windows users should run it
inside an x86-64 WSL2 Linux distribution. Native Windows, Alpine/musl Linux,
and Linux ARM are not currently supported.

For repository development, request an editable installation:

```console
bash tools/setup/install_postgwas.sh --editable
```

The installer refuses to alter an existing environment unless the user makes
that choice explicitly. Update an existing installation or use another name
with:

```console
bash tools/setup/install_postgwas.sh --update
bash tools/setup/install_postgwas.sh --name postgwas-v2-local --editable
```

For troubleshooting or controlled automation, the equivalent manual commands
are:

```console
mamba env create -f environment.yml
conda activate postgwas
python -m pip install --no-deps --no-build-isolation .
bash tools/setup/install_bcftools_liftover.sh
```

Use `conda env create` instead of `mamba env create` when Mamba is unavailable.
The environment file is the dependency source of truth; `--no-deps` prevents a
second pip resolution from replacing its tested numerical and single-cell
versions.

The environment isolates its R library from per-user R packages. This prevents
an incompatible package in `~/Library/R` or an equivalent user library from
shadowing the Conda-tested R stack.

## Verify the software layer

Run all checks in the activated environment:

```console
python --version
Rscript --version
python -m pip check
postgwas --help
postgwas config --help
bcftools --version
bcftools plugin -l | grep -x liftover
tabix --version
plink --version
scdrs --help
```

Every command must exit with status 0. `postgwas --help` proves that the public
interface imports; it does not prove that external software or scientific
resources for a selected module are ready.

For a development checkout, run the installation-focused regression tests:

```console
PYTHONDONTWRITEBYTECODE=1 \
python -m pytest -q -p no:cacheprovider \
  tests/test_cli_help.py tests/test_kpops_caldera.py
```

If scanpy reports that a Numba cache has no writable locator on a restricted
cluster, assign a writable job-local cache before starting PostGWAS:

```console
export NUMBA_CACHE_DIR="${TMPDIR:-/tmp}/postgwas-numba-${USER}"
mkdir -p "$NUMBA_CACHE_DIR"
```

## Module requirements

“In environment” below means installed by `environment.yml` and the liftover
setup command. Reference files are never inferred: build, chromosome naming,
allele convention, and population must agree with the study.

| Module or engine | Additional software beyond the local environment | Scientific resources or access |
|---|---|---|
| `harmonisation` | None after the bcftools liftover setup | Reference FASTA and index, dbSNP and frequency VCFs plus indexes, chain files for build conversion, and configured build-check files |
| `sumstat_filter`, `formatter`, `qc` | None; bcftools, tabix, bash, and Python are in the environment | Harmonised GWAS-VCF; formatter targets may need downstream metadata |
| `annot_ldblock` | None | Build- and population-matched LD-block BED files |
| `ld_clump` | None; bcftools and tabix are in the environment | Indexed population-specific pairwise-LD tables and LD blocks |
| `imputation` (`pred_ld`) | None | PRED-LD reference panel matching build and ancestry, plus the PostGWAS resource root |
| `manhattan` | None; R and plotting packages are in the environment | Harmonised GWAS-VCF |
| `finemap` (`susie`) | None; PLINK 1.9, R, and `susieR` are in the environment | PLINK genotype reference matching build, ancestry, and allele coding |
| `finemap` (`finemap`) | PLINK 2, BGENIX, LDstore 2, and FINEMAP 1.4.2 | Matching PLINK genotype reference; these executables are not native Apple Silicon packages |
| `magma`, `magmacovar`, single-cell MAGMA mode | MAGMA 1.10 | PLINK LD reference, MAGMA gene locations, and analysis-specific gene sets or covariates |
| `gcta_cojo`, `gcta_gene` | GCTA (`gcta64`); bcftools is already installed | Matching PLINK reference, gene list, and optional GMT/set files required by the selected method |
| `heritability`, single-cell LDSC mode | `ldsc.py` and `munge_sumstats.py` from one compatible LDSC installation | Matching LD scores, regression weights, HapMap3 allele list, and cell-type `.ldcts` resources when applicable |
| `pops` | None | PoPS feature matrix chunks, row/column files, and gene annotation |
| `kpops` | Pinned `k-pops.py`, installed by the repository preparation script | K-POPS kernel and gene annotation matching the MAGMA/PoPS build |
| `caldera` | Pinned CALDERA repository and R adapter | CALDERA model/data bundle, PoPS results, and fine-mapped credible sets; pipeline mode currently requires GRCh37 |
| `flames` | Network access to Ensembl VEP, or a configured local VEP/CADD path where supported | FLAMES annotations/model plus compatible fine-mapping, MAGMA, MAGMAcovar, and PoPS results |
| `mixer` | MiXeR/GSA-MiXeR runtime, normally through its configured Linux container | Matching BIM/LD patterns and GO tables for GSA |
| `single_cell` (`scdrs`) | None; scDRS and the Python stack are installed | H5AD atlas, gene-ID mapping, gene set, covariates, and configured cell/group annotations |
| `pathway_enrichment` | Dedicated enrichment Python/R environment expected by the current dispatcher | Network access, BioGRID key, DAVID-registered email, and DSigDB GMT for DSigDB |

## Install external programs

Use the executable names in `src/postgwas/config/defaults/resources.yaml`, or
set absolute paths under `resources.executables` in the run configuration.
Never point a configuration at a different version without recording and
validating that version.

- [bcftools and HTSlib](https://github.com/samtools/bcftools/releases) are
  installed by Conda. The repository setup script builds the pinned liftover
  plugin from the upstream source revision.
- [PLINK 1.9 and PLINK 2](https://www.cog-genomics.org/plink/) are separate
  programs. PLINK 1.9 is installed on supported local platforms; PLINK 2 is a
  Linux-selected Conda dependency and is additionally required by the FINEMAP
  engine.
- [GCTA](https://yanglab.westlake.edu.cn/software/gcta/),
  [MAGMA](https://cncr.nl/research/magma/), and
  [FINEMAP/LDstore](https://www.christianbenner.com/) must be obtained under
  their upstream terms. Do not redistribute a locally downloaded MAGMA binary
  without the authors' permission.
- [LDSC](https://github.com/CBIIT/ldsc) has its own tightly constrained Python
  dependencies. Prefer a separate LDSC environment and configure the absolute
  `ldsc.py` and `munge_sumstats.py` paths instead of changing the main
  PostGWAS numerical stack.
- [MiXeR](https://github.com/precimed/mixer) is configured to use a Linux
  container by default. A locally compiled runtime must supply both the MiXeR
  analysis and figure scripts plus its native library.

After each manual installation, verify the exact executable that PostGWAS will
resolve:

```console
command -v plink2 gcta64 magma finemap ldstore bgenix \
  ldsc.py munge_sumstats.py
```

Missing output means every module that depends on that program remains
unavailable.

## Prepare K-POPS and CALDERA

The preparation script downloads checksum-pinned upstream commits. It also
expects the PoPS GRCh37 gene annotation named by the script to already exist in
the resource root.

```console
conda activate postgwas
bash tools/resource_preparation/prepare_kpops_caldera_resources.sh \
  --resource-root /absolute/path/to/resources \
  --python "$CONDA_PREFIX/bin/python"
```

Add `--prepare-linear-kernel` only after the documented PoPS munged feature
chunks exist. The script validates checksums and refuses to replace an
unrecognised existing installation.

## Pathway-enrichment installation status

The current pathway-enrichment dispatcher switches to a separate micromamba
environment named `enricher`. Its complete dependency set is not installed by
`environment.yml`, and the Apple Silicon channel solve does not currently
provide `r-webgestaltr`. Therefore a clean native macOS installation is **not
ready for `pathway_enrichment`** even when `postgwas pathway_enrichment --help`
works. Use the Linux container after validating it, or complete and test the
dedicated environment before relying on this module.

## Container installation

The Dockerfile attempts to assemble the broadest Linux x86-64 tool stack,
including upstream binaries that are unavailable natively on Apple Silicon.
Build and smoke-test it from the repository root:

```console
docker build --platform linux/amd64 -t postgwas:local .
docker run --rm --platform linux/amd64 postgwas:local postgwas --help
```

Mount study data and resources rather than copying them into the image:

```console
docker run --rm --platform linux/amd64 \
  -v /absolute/path/to/project:/work \
  -v /absolute/path/to/resources:/resources:ro \
  postgwas:local \
  postgwas config validate --config /work/run.yaml
```

Use `/work` and `/resources` paths inside container configuration files. Review
the Dockerfile, upstream licences, downloaded versions, and checksums before
publishing an image. A successful `postgwas --help` smoke test still does not
validate scientific resources.

## Scientific resource readiness

No universal reference bundle ships with the package. At present the
repository documentation does not publish a complete authoritative bundle URL
and checksum manifest, so a new user cannot reproduce every module from the
checkout alone. Follow [Resource Setup](resource-setup.md) and
[Reference Resources](../reference/reference-resources.md), and retain for each
file:

- source URL and retrieval date;
- release/version and checksum;
- genome build and chromosome naming;
- population/reference-panel provenance;
- allele and variant-ID convention; and
- every preparation command.

Do not substitute a reference from a different build or ancestry simply to
make preflight pass. That can change LD, allele alignment, fine-mapping,
heritability, gene association, and prioritisation results.

## Installation acceptance checklist

A user installation is ready for a selected analysis only when all of these are
true:

- the clean environment solve and package installation complete without error;
- `python -m pip check` reports no broken requirements;
- the public command and selected module help pages exit successfully;
- every required executable resolves to the intended version;
- the bcftools liftover plugin is present when build conversion is requested;
- every required R/Python library imports in the activated environment;
- all input and reference files exist with their indexes;
- build, population, alleles, IDs, and sample-size definitions are compatible;
- the module preflight succeeds; and
- a small representative run produces and validates every required output.

Configuration validation alone is not an end-to-end installation test. Never
interpret a result until the module's required outputs, QC, log, and completion
status all confirm success.
