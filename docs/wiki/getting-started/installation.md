# Installation

PostGWAS has three separate readiness layers:

1. the PostGWAS package and its Python/R libraries;
2. third-party command-line programs used by selected modules; and
3. genome-build-, population-, and method-specific reference resources.

The complete Mamba installer installs and checks the first two layers. A
successful installation is not a test of every analysis on your study data.
Reference resources remain explicit inputs because build, ancestry, release,
allele convention, credentials, and licensing differ by analysis. A module is
ready only after every item in its row of the
[module requirements table](#module-requirements) is present and compatible.

## Supported installation paths

The table describes local Mamba/Conda installations. The repository also supplies
a [Linux-amd64 Docker build recipe](#docker-installation), with separate
platform requirements and validation limits described below.

| Host | Portable environment | Complete Mamba environment |
|---|---|---|
| glibc Linux x86-64 | Supported | Supported and recommended |
| Linux ARM64 | Not supported | Not supported |
| macOS Apple Silicon | Supported for the portable stack | Supported with Rosetta 2 and Apple Command Line Tools |
| macOS Intel | Supported for the portable stack | Not supported |
| Windows | Not supported natively | Use a glibc Linux x86-64 WSL2 distribution |

The environment uses [CEP-24 platform selectors](https://conda.org/learn/ceps/cep-0024/).
Use a current Mamba or Conda release that understands dictionary selectors;
older clients may reject entries such as `sel(linux)`.

## Prerequisites for local installation

- Git and a complete checkout of this repository.
- [Miniforge](https://github.com/conda-forge/miniforge) with Mamba for the
  complete installation.
- Apple-silicon macOS with Apple Command Line Tools and Rosetta 2, or a
  glibc-based Linux x86-64 distribution such as Ubuntu, Debian, Rocky Linux,
  or Fedora, for `--all-tools`.
- Enough storage for the main environment, isolated compatibility runtimes, source
  builds, and temporary downloads. Scientific reference data require
  substantially more space.
- Outbound HTTPS access. FINEMAP and LDstore are served by their official site
  over HTTP; the pinned SHA-256 checksums remain mandatory for those archives.

Do not use the system Python. The all-module environment is pinned to Python
3.10 because the current scanpy/scDRS/Numba stack is not compatible with
Python 3.12 or later.

## Install the complete software stack

After installing and initialising Miniforge with Mamba, clone the repository
and run the installer on Apple-silicon macOS or glibc Linux x86-64:

```console
git clone https://github.com/JIBINJOHNV/postgwas_v2.git
cd postgwas_v2
bash tools/setup/install_postgwas.sh \
  --all-tools \
  --name postgwas
conda activate postgwas
```

`--all-tools` requires Mamba and always creates a new environment. It refuses
to alter an existing environment. Choose another name if `postgwas` already
exists:

```console
bash tools/setup/install_postgwas.sh \
  --all-tools \
  --name postgwas-study
```

Source compilation uses two jobs by default. Set an explicit bounded value for
the host or scheduler allocation:

```console
bash tools/setup/install_postgwas.sh \
  --all-tools \
  --name postgwas-study \
  --jobs 8
```

The installer performs these operations in order:

1. creates the named Mamba environment from `environment.yml`;
2. installs the checkout without a second dependency resolution;
3. downloads, checksum-verifies, builds, and installs bcftools plus its
   liftover plugin;
4. downloads verified MAGMA, FINEMAP, and LDstore distributions and, on
   macOS, official PLINK 2 and GCTA distributions plus verified BGENIX source;
5. installs immutable LDSC, K-POPS, CALDERA, and MiXeR upstream revisions;
6. creates dependency-isolated LDSC, pathway-enrichment, and VEP runtimes
   beneath the main environment, plus Intel-tool and MiXeR runtimes on macOS;
7. exposes their commands through the activated PostGWAS environment; and
8. runs the complete installation verifier.

The installation exits non-zero when a download, checksum, environment solve,
source build, import, executable, native library, or version check fails. The
version, URL, commit, and checksum source of truth is
`tools/setup/software_versions.env`.

FLaMEs 1.1.2-compatible code, its model, and PoPS are packaged PostGWAS
modules, so they do not require a second external executable. The verifier
imports both modules and checks the surrounding numerical stack. Analysis
still requires the method-specific data listed below.

No separate `gwas2vcf` package or executable is downloaded. The protected
adapter code already shipped inside PostGWAS remains available to the
harmonisation implementation.

## Why some runtimes are isolated

Users activate only the main environment and run ordinary commands such as
`postgwas heritability` or `postgwas pathway_enrichment`. Internally, the
installer places these runtimes under `$CONDA_PREFIX/share/postgwas/environments`:

- LDSC uses its upstream-pinned NumPy, pandas, SciPy, pysam, and bitarray
  versions. `ldsc.py` and `munge_sumstats.py` in the main environment point to
  that compatible runtime.
- Pathway enrichment uses an isolated Python/R environment containing rpy2,
  OmniPath, WebGestaltR, and the remaining provider clients. The dispatcher
  invokes its installed Python directly.
- Ensembl VEP uses a separate Bioconda environment. The `vep` command is
  exposed in the main environment; VEP cache and FASTA data are not installed.
- On Apple silicon, FINEMAP and LDstore use a small Intel library runtime, and
  MiXeR uses its pinned Intel Python/native stack. Their public commands remain
  in the main environment and dispatch through Rosetta automatically.

Combining these dependency sets with the main single-cell and numerical stack
would make the environment inconsistent. This is one user-facing activation,
not one dependency environment. Isolation is an installation detail, not a
second user activation step.

Each created environment clears inherited `PYTHONPATH` and disables per-user
Python site packages. The installation verifier confirms that PostGWAS,
FLaMEs, and PoPS resolve from the new environment or the selected current
checkout, never from a sibling or retired repository. Old checkouts may remain
on disk, but the installed commands do not search or execute them.

## Portable installation

Use portable mode when the complete third-party stack is not needed:

```console
bash tools/setup/install_postgwas.sh --name postgwas-portable
conda activate postgwas-portable
```

The script uses Mamba when available and otherwise Conda. Portable mode
supports macOS arm64/x86-64 and glibc Linux x86-64. It installs the package,
Python/R dependencies, portable Conda programs, bcftools, and the liftover
plugin. It does not install MAGMA, LDSC, FINEMAP/LDstore, MiXeR, K-POPS,
CALDERA, VEP, or pathway enrichment.

For repository development, request an editable installation:

```console
bash tools/setup/install_postgwas.sh \
  --name postgwas-development \
  --editable
```

Portable mode alone may update its own existing environment when explicitly
requested:

```console
bash tools/setup/install_postgwas.sh \
  --name postgwas-portable \
  --update
```

The equivalent portable manual commands are:

```console
mamba env create --name postgwas-portable --file environment.yml
conda activate postgwas-portable
python -m pip install --no-deps --no-build-isolation .
bash tools/setup/install_bcftools_liftover.sh
```

The environment isolates its R library from per-user R packages, preventing an
incompatible package in `~/Library/R` or an equivalent user library from
shadowing the tested R stack.

## Docker installation

The repository includes a [Dockerfile](https://github.com/JIBINJOHNV/postgwas_v2/blob/main/Dockerfile)
that builds PostGWAS, its external programs and isolated runtimes into one local
image. You do not need Mamba, Python or R installed on the host for this route.
Reference panels, annotation/cache data and service credentials are not included.

### Prepare Docker and choose the platform

Install and start [Docker Desktop on macOS](https://docs.docker.com/desktop/setup/install/mac-install/)
or [Docker Engine on Linux](https://docs.docker.com/engine/install/), then check
that `docker info` succeeds. Git is needed to obtain the build sources.

The recipe targets **Linux x86-64 (`linux/amd64`)**, including its MiXeR
compiler settings and FINEMAP/LDstore binaries. It is not a native ARM64 build.
Apple-silicon Macs require amd64 emulation; Docker supports this, but compilation
and analyses can be slower and need platform-specific testing. See
[Docker's emulation guidance](https://docs.docker.com/build/building/multi-platform/).
Allocate sufficient RAM, CPU and disk space to Docker for the build and study;
PostGWAS compute settings must fit within those allocations.

### Build and check the image

Use a clean checkout with no private study files, credentials or local notes.
The Dockerfile copies the build context into the image: `.gitignore` does not
replace [Docker's `.dockerignore` rules](https://docs.docker.com/build/concepts/context/#dockerignore-files).
Keep study data outside the checkout and mount it only when running analyses.

```console
git clone https://github.com/JIBINJOHNV/postgwas_v2.git
cd postgwas_v2
docker build --platform linux/amd64 --load --tag postgwas:local .
docker run --rm --platform linux/amd64 postgwas:local postgwas --help
```

`postgwas:local` is the image you build, not a published image to pull. Record
the source commit and resulting image ID with the analysis. Rebuild after
updating the checkout; an existing image does not change when local files change.
The Dockerfile has no PostGWAS entrypoint, so commands after the image name must
include `postgwas`, as above.

The build runs `tools/setup/verify_all_tools.sh`, checking program availability,
selected imports, dependency consistency, versions and required model/library
files. You can repeat those checks after a successful build:

```console
docker run --rm --platform linux/amd64 postgwas:local \
  bash /opt/postgwas/tools/setup/verify_all_tools.sh
```

Current installation CI tests the Mamba route, not Docker builds or full module
analyses in the image. A Docker build recipe and passing software checks are
not a guarantee that every analysis works on every host; validate a representative
study before a full run. Review all third-party licences before sharing the
image, especially the [MAGMA redistribution terms](https://cncr.nl/research/magma/).

### Mount files and run an analysis

Replace the three host paths below with existing directories. The input folder
must contain the indexed PostGWAS-harmonised VCF used in the example. Prepare
reference files and their indexes before mounting them read-only. Create a
separate writable output directory on the host first.

Start an interactive container shell from a macOS or Linux terminal:

```console
docker run --rm -it --platform linux/amd64 \
  --user "$(id -u):$(id -g)" \
  --mount "type=bind,src=/absolute/path/to/inputs,dst=/input,readonly" \
  --mount "type=bind,src=/absolute/path/to/references,dst=/reference,readonly" \
  --mount "type=bind,src=/absolute/path/to/outputs,dst=/output" \
  --workdir /output \
  postgwas:local bash
```

Inside that shell, run ordinary PostGWAS commands. For example, assess the VCF
through pipeline mode (QC itself does not use the mounted reference folder):

```console
postgwas pipeline --modules qc_summary \
  --vcf /input/STUDY_GRCh37_merged.vcf.gz \
  --dataset-id STUDY \
  --output-directory /output/qc
```

Use `/input`, `/reference` and `/output` paths in commands, sample sheets and
run YAML. Host paths such as `/Users/...` are unavailable unless mounted at
those same locations. Use absolute container paths to avoid differences in
command, sample-sheet and configuration path resolution.
The host user/group selection avoids root-owned output on Linux; ensure that
user can read the inputs and write the output folder. Docker Desktop must be
allowed to share the selected folders. See [Docker bind mounts](https://docs.docker.com/engine/storage/bind-mounts/).

Results and logs written under `/output` remain in the host output directory
when you exit. `--rm` removes the container, so do not store required results
only in its unmounted filesystem. Use the module guides for additional inputs,
references and parameters; the image does not remove those requirements.

### Full image versus MiXeR's Docker backend

The full PostGWAS image already includes a MiXeR runtime. Use
`--mixer-backend native` when selecting MiXeR inside this image; “native” here
means inside the running container. Do not request nested Docker execution.

A host-installed PostGWAS can separately select `--mixer-backend docker` to run
only MiXeR in the configured GSA-MiXeR image. That backend does not install or
containerise the other PostGWAS modules. See the [MiXeR guide](../modules/mixer.md).

## Verify the local installation

The complete installer runs this automatically after every tool is installed:

```console
bash tools/setup/verify_mamba_installation.sh
```

Useful manual checks are:

```console
python -m pip check
postgwas --help
postgwas config --help
bcftools --version
bcftools plugin -l | grep -x liftover
magma --version
ldsc.py --help
munge_sumstats.py --help
command -v finemap ldstore bgenix
k-pops.py --help
mixer.py --version
vep --help
```

Every command must resolve inside the activated environment and exit
successfully. `postgwas --help` only proves that the public interface imports;
it does not prove that scientific resources for a selected module are ready.

The test suite is retained locally in maintainer copies and is not included in
a public Git clone. If your development checkout has that local-only suite,
you can also run the focused installation contracts below. Users without it
should use the shipped verification scripts and manual checks above; a missing
`tests/` directory is not an installation failure.

```console
PYTHONDONTWRITEBYTECODE=1 \
python -m pytest -q -p no:cacheprovider tests/test_installation_script.py
```

If scanpy reports that a Numba cache has no writable locator on a restricted
cluster, assign a writable job-local cache before starting PostGWAS:

```console
export NUMBA_CACHE_DIR="${TMPDIR:-/tmp}/postgwas-numba-${USER}"
mkdir -p "$NUMBA_CACHE_DIR"
```

## Module requirements

“Installed” below means supplied by the complete Mamba installer. Reference
files are never inferred: build, chromosome naming, allele convention, and
population must agree with the study.

| Module or engine | Software status after `--all-tools` | Scientific resources or access |
|---|---|---|
| `harmonisation` | bcftools, tabix, and the liftover plugin installed | Reference FASTA/index, dbSNP and frequency resources/indexes, GFF annotation, chain files, and configured build-check files |
| `sumstat_filter`, `formatter`, `qc` | Installed | Harmonised GWAS-VCF; formatter targets may require downstream metadata |
| `annot_ldblock` | Installed | Build- and population-matched LD-block BED files |
| `ld_clump` (`standard`) | Installed | Prepared manifest-validated, indexed pairwise-LD reference and variant inventories |
| `ld_clump` (`region`) | Installed | Population-specific LD-block annotation in the VCF, or BED files for pipeline preparation |
| `ld_clump` (`cojo-slct`) | GCTA installed | Matching PLINK BED/BIM/FAM reference; see the method-specific clumping guide |
| `imputation` (`pred_ld`) | Installed | Matching PRED-LD reference plus the full harmonisation resource tree, because imputed output is re-harmonised |
| `manhattan` | R and plotting packages installed | Harmonised GWAS-VCF |
| `finemap` (`susie`) | PLINK 1.9, R, and `susieR` installed | Matching PLINK genotypes; pipeline mode also needs the standard-clumping pairwise-LD bundle |
| `finemap` (`finemap`) | PLINK 2, BGENIX, LDstore 2, and FINEMAP 1.4.2 installed | Matching PLINK genotypes; pipeline mode also needs the standard-clumping pairwise-LD bundle |
| `magma`, `magmacovar`, single-cell MAGMA | MAGMA 1.10 installed | PLINK LD reference, gene locations, and analysis-specific gene sets/covariates |
| `gcta_cojo` | GCTA installed | Matching PLINK reference and method-specific SNP lists when required |
| `gcta_gene` | GCTA installed | Matching PLINK reference; gene coordinates for gene tests or GMT conversion, a native set list or GMT for set tests, and no gene list for fixed-segment tests |
| `heritability`, single-cell LDSC | Isolated CBIIT LDSC installed | Matching LD scores, regression weights, HapMap3 alleles, and optional `.ldcts` resources |
| `pops` | Packaged module and dependencies installed | PoPS feature chunks, row/column files, and gene annotation |
| `kpops` | Pinned `k-pops.py` installed | K-POPS kernel and gene annotation matching the MAGMA/PoPS build |
| `caldera` | Pinned CALDERA runner/model installed | PoPS results and fine-mapped credible sets; pipeline mode currently requires GRCh37 |
| `flames` | Packaged code/model and VEP command installed; API mode remains available | FLAMES annotation bundle plus compatible fine-mapping, MAGMA, MAGMAcovar, and PoPS results; local modes require VEP/CADD data |
| `mixer` | Native MiXeR/GSA-MiXeR runtime installed | Matching BIM/LD patterns and GO tables for GSA |
| `single_cell` (`scdrs`) | scDRS and Python stack installed | H5AD atlas, gene-ID mapping, gene set, covariates, and configured annotations |
| `pathway_enrichment` | Isolated Python/R provider runtime installed | Network access, BioGRID key, DAVID-registered email, and DSigDB GMT |

## External program provenance and licensing

- [bcftools and HTSlib](https://github.com/samtools/bcftools/releases) are
  built from the pinned official release. The liftover plugin is built from a
  pinned `freeseek/score` revision.
- [PLINK 1.9 and PLINK 2](https://www.cog-genomics.org/plink/) are separate
  programs installed through the canonical Mamba environment.
- [GCTA](https://yanglab.westlake.edu.cn/software/gcta/) and BGENIX use their
  Bioconda packages on Linux; on macOS, GCTA uses its official arm64 archive
  and BGENIX is checksum-verified and built from its official source release.
- [MAGMA](https://cncr.nl/research/magma/) 1.10 is downloaded from the official
  CNCR link. MAGMA states that post-v1.0 binaries and source may not be
  redistributed or modified, so the repository stores only the downloader,
  version, and checksum—not the binary.
- [LDSC](https://github.com/CBIIT/ldsc) is installed from an immutable commit in
  its own environment because its numerical requirements conflict with the
  main stack.
- [FINEMAP and LDstore](https://www.christianbenner.com/) are downloaded from
  their official distributions and checksum-verified.
- [FLaMEs](https://github.com/Marijn-Schipper/FLAMES) and
  [PoPS](https://github.com/FinucaneLab/pops) are integrated package modules.
- [MiXeR](https://github.com/precimed/gsa-mixer),
  [K-POPS](https://github.com/JasonTan-code/k-pops), and
  [CALDERA](https://github.com/kheilbron/caldera) are installed from immutable
  upstream revisions.

Review every upstream licence before copying or redistributing an environment.

## Prepare K-POPS scientific resources

The complete installer supplies the K-POPS command and CALDERA software. To
validate and prepare K-POPS scientific resources later, provide the resource
root containing the expected PoPS GRCh37 gene annotation:

```console
conda activate postgwas
bash tools/resource_preparation/prepare_kpops_caldera_resources.sh \
  --resource-root /absolute/path/to/resources \
  --python "$CONDA_PREFIX/bin/python"
```

Add `--prepare-linear-kernel` only after the documented PoPS munged feature
chunks exist. The script validates checksums and refuses to replace an
unrecognised existing installation.

## Scientific resource readiness

No universal reference bundle ships with the package. The repository does not
currently publish one complete authoritative bundle URL and checksum manifest,
so installing software alone cannot reproduce every module. Follow
[Resource Setup](resource-setup.md) and
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

A user installation is ready for a selected analysis only when all of these
are true:

- the clean environment solve and software installation complete without error;
- `python -m pip check` reports no broken requirements in every relevant runtime;
- the public command and selected module help pages exit successfully;
- every required executable resolves to the intended version;
- the bcftools liftover plugin is present when build conversion is requested;
- every required R/Python library imports;
- all input and reference files exist with their indexes;
- build, population, alleles, IDs, and sample-size definitions are compatible;
- the module preflight succeeds; and
- a representative run produces and validates every required output.

Configuration validation alone is not an end-to-end installation test. Never
interpret a result until the module's required outputs, QC, log, and completion
status all confirm success.

## Next steps

1. [Prepare the references required by your selected analyses](resource-setup.md).
2. [Follow the connected harmonisation, QC, and MAGMA tutorial](quick-start.md).
3. If a check fails, use [Troubleshooting](../help/troubleshooting.md) before
   starting a full analysis.
