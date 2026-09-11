# =====================================================================
# MiXeR build stage
# =====================================================================
# This temporary stage compiles MiXeR. Boost, source code, Git history,
# tests and build files are not copied into the final PostGWAS image.
FROM mambaorg/micromamba:2.8.1@sha256:fb18405d6004af757a38ec498a078240b4fd5549146990a484c28bb7e78aace4 AS mixer-builder

USER root

ARG GSA_MIXER_COMMIT=cf65c57d5d1ad76597db1d4fa3907f1d711d39e7
ARG MIXER_BUILD_JOBS=4

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential ca-certificates cmake git wget && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /tmp

RUN wget -q \
    https://precimed.s3.eu-west-1.amazonaws.com/gsa-mixer/boost_1_69_0.tar.gz && \
    tar -xzf boost_1_69_0.tar.gz && \
    cd boost_1_69_0 && \
    ./bootstrap.sh \
    --with-libraries=program_options,filesystem,system,date_time && \
    ./b2 \
    -j"${MIXER_BUILD_JOBS}" \
    cxxflags=-fPIC \
    variant=release \
    link=static \
    threading=multi \
    --with-program_options \
    --with-filesystem \
    --with-system \
    --with-date_time

# Replace host-specific CPU compilation with the minimum SIMD features
# required by MiXeR's FastDifferentialCoding implementation. The pinned
# TurboPFor source relies on removed C implicit-int behaviour in two decoder
# macros, and its bundled NLopt source omits the standard allocation header.
# These declarations preserve the original arithmetic and control flow while
# allowing modern GCC to compile the pinned source.
RUN mkdir -p /tmp/gsa-mixer && \
    git -C /tmp/gsa-mixer init && \
    git -C /tmp/gsa-mixer remote add origin \
    https://github.com/precimed/gsa-mixer.git && \
    git -C /tmp/gsa-mixer fetch --depth 1 origin "${GSA_MIXER_COMMIT}" && \
    git -C /tmp/gsa-mixer checkout --detach FETCH_HEAD && \
    test "$(git -C /tmp/gsa-mixer rev-parse HEAD)" = "${GSA_MIXER_COMMIT}" && \
    rm -rf /tmp/gsa-mixer/.git && \
    sed -i 's/-march=native/-mssse3 -msse4.1/g' \
    /tmp/gsa-mixer/src/CMakeLists.txt && \
    mixer_bitutil=/tmp/gsa-mixer/src/TurboPFor/bitutil.c && \
    test "$(grep -c 'const _md = _md_;' "${mixer_bitutil}")" -eq 2 && \
    sed -i 's/const _md = _md_;/const _t_ _md = _md_;/g' \
    "${mixer_bitutil}" && \
    test "$(grep -c 'const _t_ _md = _md_;' "${mixer_bitutil}")" -eq 2 && \
    mixer_nlopt_stop=/tmp/gsa-mixer/src/nlopt/stop.c && \
    test "$(grep -c '^#include <stdarg.h>$' "${mixer_nlopt_stop}")" -eq 1 && \
    test "$(grep -c '^#include <stdlib.h>$' "${mixer_nlopt_stop}")" -eq 0 && \
    sed -i '/^#include <stdarg.h>$/a #include <stdlib.h>' \
    "${mixer_nlopt_stop}" && \
    test "$(grep -c '^#include <stdlib.h>$' "${mixer_nlopt_stop}")" -eq 1 && \
    cmake \
    -S /tmp/gsa-mixer/src \
    -B /tmp/gsa-mixer/src/build \
    -DBoost_NO_BOOST_CMAKE=ON \
    -DBOOST_ROOT=/tmp/boost_1_69_0 \
    -DBoost_USE_STATIC_LIBS=ON && \
    cmake \
    --build /tmp/gsa-mixer/src/build \
    --target bgmg \
    --parallel "${MIXER_BUILD_JOBS}"

# Copy only the files needed to run MiXeR.
RUN mkdir -p /mixer-runtime/precimed /mixer-runtime/lib && \
    cp -a /tmp/gsa-mixer/precimed/. /mixer-runtime/precimed/ && \
    cp /tmp/gsa-mixer/src/build/lib/libbgmg.so /mixer-runtime/lib/ && \
    cp /tmp/gsa-mixer/LICENSE /mixer-runtime/LICENSE && \
    rm -rf \
    /mixer-runtime/precimed/mixer-test \
    /mixer-runtime/precimed/__pycache__ && \
    find /mixer-runtime/precimed -type f \
    \( -name '*.pyc' -o -name '*.pyo' \) -delete && \
    strip --strip-unneeded /mixer-runtime/lib/libbgmg.so


# =====================================================================
# Base Image
# =====================================================================
FROM mambaorg/micromamba:2.8.1@sha256:fb18405d6004af757a38ec498a078240b4fd5549146990a484c28bb7e78aace4

USER root

# =====================================================================
# System dependencies
# =====================================================================
# We use apt-get for lightweight system libraries.
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    gcc g++ make wget git curl bzip2 unzip zip \
    libcurl4-openssl-dev libsuitesparse-dev \
    zlib1g-dev libbz2-dev liblzma-dev \
    ca-certificates pkg-config \
    libssl-dev libxml2-dev pigz \
    libharfbuzz-dev \
    gnupg \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff-dev \
    libjpeg-dev \
    libfontconfig1-dev \
    libgit2-dev \
    libeigen3-dev \
    procps && \
    rm -rf /var/lib/apt/lists/*

# =====================================================================
# 1. MAIN ENV: PostGWAS
# =====================================================================
COPY environment.yml /tmp/environment.yml
COPY tools/setup/install_bcftools_liftover.sh /tmp/install_bcftools_liftover.sh
COPY tools/setup/software_versions.env /tmp/software_versions.env
COPY tools/setup/container/ldsc_wrapper.sh /usr/local/libexec/postgwas/ldsc_wrapper.sh

# Create env and clean immediately.
RUN micromamba create -y -n postgwas -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Add Linux-selected tools plus the R and MiXeR runtime dependencies.
RUN linux_packages="$( \
    awk '$1 == "-" && $2 == "sel(linux):" { print $3 }' \
    /tmp/environment.yml \
    )" && \
    test -n "${linux_packages}" && \
    micromamba install -y -n postgwas -c conda-forge -c bioconda \
    cmake \
    mkl \
    mkl-include \
    sysroot_linux-64 \
    gcc_linux-64 \
    gxx_linux-64 \
    "r-base=4.3.*" \
    r-xml2 \
    libxml2 \
    r-matrix \
    r-data.table \
    r-survey \
    intervaltree \
    "matplotlib<3.10" \
    matplotlib-venn \
    numdifftools \
    "pandas<3" \
    six \
    ${linux_packages} && \
    micromamba clean --all --yes

ENV PATH="/opt/conda/envs/postgwas/bin:$PATH"

RUN micromamba run -n postgwas \
    bash /tmp/install_bcftools_liftover.sh && \
    rm /tmp/install_bcftools_liftover.sh && \
    micromamba clean --all --yes

ENV BCFTOOLS_PLUGINS="/opt/conda/envs/postgwas/libexec/bcftools"

# =====================================================================
# 2. VEP
# =====================================================================
RUN micromamba create -y -n vep -c conda-forge -c bioconda ensembl-vep=113 && \
    micromamba clean --all --yes

# LDSC pins older numerical packages that conflict with the main PostGWAS
# environment. Keep it isolated, and expose stable command wrappers on PATH.
ARG LDSC_REPOSITORY=https://github.com/CBIIT/ldsc.git
ARG LDSC_COMMIT=6c673952cee74bd5c57aef1555a03b1c015399a0
RUN micromamba create -y -n ldsc -c conda-forge -c bioconda \
    python=3.10 \
    pip \
    bitarray=2.6.0 \
    nose=1.3.7 \
    numpy=1.23.3 \
    pandas=1.5.0 \
    pybedtools=0.9.1 \
    pysam=0.19.1 \
    python-dateutil=2.8.2 \
    pytz=2022.4 \
    scipy=1.9.2 \
    six=1.16.0 && \
    mkdir -p /opt/ldsc && \
    git -C /opt/ldsc init && \
    git -C /opt/ldsc remote add origin "${LDSC_REPOSITORY}" && \
    git -C /opt/ldsc fetch --depth 1 origin "${LDSC_COMMIT}" && \
    git -C /opt/ldsc checkout --detach FETCH_HEAD && \
    test "$(git -C /opt/ldsc rev-parse HEAD)" = "${LDSC_COMMIT}" && \
    rm -rf /opt/ldsc/.git && \
    micromamba run -n ldsc \
    python -m pip install --no-deps --no-cache-dir /opt/ldsc && \
    micromamba run -n ldsc python -m pip check && \
    micromamba run -n ldsc ldsc.py --help >/dev/null && \
    micromamba run -n ldsc munge_sumstats.py --help >/dev/null && \
    chmod 0755 /usr/local/libexec/postgwas/ldsc_wrapper.sh && \
    ln -s /usr/local/libexec/postgwas/ldsc_wrapper.sh /usr/local/bin/ldsc.py && \
    ln -s /usr/local/libexec/postgwas/ldsc_wrapper.sh /usr/local/bin/munge_sumstats.py && \
    micromamba clean --all --yes

# K-POPS and CALDERA run in the same PostGWAS conda environment. Both
# repositories are pinned to immutable commits because neither project
# currently publishes a versioned conda package.
ARG KPOPS_REPOSITORY=https://github.com/JasonTan-code/k-pops.git
ARG KPOPS_COMMIT=8acd49ed8c96565b17c2997420f514f24b65e096
ARG CALDERA_REPOSITORY=https://github.com/kheilbron/caldera.git
ARG CALDERA_COMMIT=81a8a0308741ae986660f711bbf6a7abbd3bca19
RUN mkdir -p /opt/kpops /opt/caldera && \
    git -C /opt/kpops init && \
    git -C /opt/kpops remote add origin "${KPOPS_REPOSITORY}" && \
    git -C /opt/kpops fetch --depth 1 origin "${KPOPS_COMMIT}" && \
    git -C /opt/kpops checkout --detach FETCH_HEAD && \
    test "$(git -C /opt/kpops rev-parse HEAD)" = "${KPOPS_COMMIT}" && \
    git -C /opt/caldera init && \
    git -C /opt/caldera remote add origin "${CALDERA_REPOSITORY}" && \
    git -C /opt/caldera fetch --depth 1 origin "${CALDERA_COMMIT}" && \
    git -C /opt/caldera checkout --detach FETCH_HEAD && \
    test "$(git -C /opt/caldera rev-parse HEAD)" = "${CALDERA_COMMIT}" && \
    rm -rf /opt/kpops/.git /opt/caldera/.git && \
    test -s /opt/kpops/k-pops.py && \
    test -s /opt/kpops/prepare_kernel.py && \
    install -m 0755 /opt/kpops/k-pops.py /opt/conda/envs/postgwas/bin/k-pops.py && \
    mkdir -p /opt/conda/envs/postgwas/share/postgwas && \
    mv /opt/caldera /opt/conda/envs/postgwas/share/postgwas/caldera && \
    test -s /opt/conda/envs/postgwas/share/postgwas/caldera/z_caldera.R && \
    test -s /opt/conda/envs/postgwas/share/postgwas/caldera/trained_models/caldera_model_no_covs.rds && \
    micromamba run -n postgwas python /opt/conda/envs/postgwas/bin/k-pops.py --help >/dev/null && \
    micromamba run -n postgwas Rscript -e \
    "parse(file='/opt/conda/envs/postgwas/share/postgwas/caldera/z_caldera.R'); library(data.table); library(dplyr)" \
    >/dev/null && \
    micromamba clean --all --yes

# =====================================================================
# 3. PATHWAY-ENRICHMENT ENVIRONMENT
# =====================================================================
# rpy2 resolves a compatible R runtime in this isolated environment.
RUN micromamba create -y -n enricher -c conda-forge -c bioconda \
    omnipath \
    gseapy \
    bgenix \
    zeep \
    rpy2 \
    pandas \
    requests \
    networkx \
    matplotlib \
    pydantic \
    pyyaml \
    rich \
    rich-argparse \
    scipy \
    "numpy<2.0" \
    r-webgestaltr \
    zip \
    gprofiler-official && \
    micromamba clean --all --yes

# =====================================================================
# Install LDStore, FINEMAP, and MAGMA
# =====================================================================
WORKDIR /tmp/install_tools

ARG LDSTORE_URL=http://www.christianbenner.com/ldstore_v2.0_x86_64.tgz
ARG LDSTORE_SHA256=be80818e2cdb5d223a15805b1fcf7569ba2b045a29bc02ca3e08af2eba95325f
ARG FINEMAP_URL=http://www.christianbenner.com/finemap_v1.4.2_x86_64.tgz
ARG FINEMAP_SHA256=3b1fc6eb3c2ccafd647b32e02d0244495cd0ade9ed7d474606c31ebf6e98b0c9
RUN curl --fail --location --silent --show-error \
    "${LDSTORE_URL}" -o ldstore_v2.0_x86_64.tgz && \
    echo "${LDSTORE_SHA256}  ldstore_v2.0_x86_64.tgz" | sha256sum --check - && \
    tar -xzf ldstore_v2.0_x86_64.tgz && \
    install -m 0755 \
    ldstore_v2.0_x86_64/ldstore_v2.0_x86_64 /usr/local/bin/ldstore && \
    curl --fail --location --silent --show-error \
    "${FINEMAP_URL}" -o finemap_v1.4.2_x86_64.tgz && \
    echo "${FINEMAP_SHA256}  finemap_v1.4.2_x86_64.tgz" | sha256sum --check - && \
    tar -xzf finemap_v1.4.2_x86_64.tgz && \
    install -m 0755 \
    finemap_v1.4.2_x86_64/finemap_v1.4.2_x86_64 /usr/local/bin/finemap && \
    rm -rf /tmp/install_tools

WORKDIR /tmp

# MAGMA v1.10 is downloaded from the static Linux link published by CNCR.
# CNCR states that post-v1.0 MAGMA binaries may not be redistributed. This
# installation is therefore appropriate for local image builds; permission
# from the MAGMA authors is required before publishing the resulting image.
ARG MAGMA_VERSION=1.10
ARG MAGMA_URL=https://vu.data.surf.nl/index.php/s/lxDgt2dNdNr6DYt/download
ARG MAGMA_SHA256=d8b20778b773f47b4fb0f1020baefebc92b11ee65946e708a618d494d6819e39
RUN mkdir -p /tmp/magma-install && \
    curl --fail --location --silent --show-error \
    "${MAGMA_URL}" -o /tmp/magma-install/magma.zip && \
    echo "${MAGMA_SHA256}  /tmp/magma-install/magma.zip" | sha256sum --check - && \
    unzip -q /tmp/magma-install/magma.zip -d /tmp/magma-install && \
    install -m 0755 /tmp/magma-install/magma /usr/local/bin/magma && \
    /usr/local/bin/magma --version > /tmp/magma-install/version.txt 2>&1
RUN grep -Eq "MAGMA version: v${MAGMA_VERSION}" /tmp/magma-install/version.txt && \
    rm -rf /tmp/magma-install

# =====================================================================
# Install MiXeR
# =====================================================================
# Only the stripped native library, Python runtime modules and licence are
# copied. MiXeR source/build folders, tests, Git data and Boost are discarded
# with the temporary mixer-builder stage.
COPY --from=mixer-builder /mixer-runtime /tools/mixer

ENV BGMG_SHARED_LIBRARY="/tools/mixer/lib/libbgmg.so"
ENV MIXER_HOME="/tools/mixer"
ENV MIXER_PY="/tools/mixer/precimed/mixer.py"
ENV MIXER_DEV_PY="/tools/mixer/precimed/mixer_dev.py"
ENV MIXER_FIGURES_PY="/tools/mixer/precimed/mixer_figures.py"
ENV POSTGWAS_ENRICHMENT_PYTHON="/opt/conda/envs/enricher/bin/python"
ENV PYTHONPATH="/tools/mixer/precimed"

RUN chmod +x \
    /tools/mixer/precimed/mixer.py \
    /tools/mixer/precimed/mixer_dev.py \
    /tools/mixer/precimed/mixer_figures.py && \
    if ldd /tools/mixer/lib/libbgmg.so | grep -q "not found"; then \
    echo "ERROR: MiXeR native library has unresolved dependencies" >&2; \
    ldd /tools/mixer/lib/libbgmg.so >&2; \
    exit 1; \
    fi && \
    micromamba run -n postgwas \
    python /tools/mixer/precimed/mixer.py --version && \
    micromamba run -n postgwas \
    python /tools/mixer/precimed/mixer.py --help >/dev/null && \
    micromamba run -n postgwas \
    python /tools/mixer/precimed/mixer_dev.py --help >/dev/null && \
    micromamba run -n postgwas \
    python /tools/mixer/precimed/mixer_figures.py --help >/dev/null

# =====================================================================
# Install PostGWAS into the main and enrichment environments
# =====================================================================
WORKDIR /opt/postgwas
COPY . /opt/postgwas

RUN micromamba run -n postgwas \
    python -m pip install --no-deps --no-build-isolation --no-cache-dir . && \
    micromamba run -n enricher \
    python -m pip install --no-deps --no-build-isolation --no-cache-dir . && \
    install -m 0755 tools/setup/mixer_wrapper.py \
    /opt/conda/envs/postgwas/bin/mixer.py && \
    install -m 0755 tools/setup/mixer_wrapper.py \
    /opt/conda/envs/postgwas/bin/mixer_dev.py && \
    install -m 0755 tools/setup/mixer_wrapper.py \
    /opt/conda/envs/postgwas/bin/mixer_figures.py && \
    micromamba run -n postgwas Rscript -e \
    "parse(file='/opt/postgwas/src/postgwas/modules/caldera/run_caldera.R')" \
    >/dev/null && \
    micromamba remove -n postgwas -y \
    cmake \
    mkl-include && \
    micromamba clean --all --yes

RUN bash /opt/postgwas/tools/setup/verify_all_tools.sh

# =====================================================================
# Final Config
# =====================================================================
ENV HOME=/tmp
RUN chown -R mambauser:mambauser /opt/postgwas
ENV PATH="/opt/conda/envs/postgwas/bin:$PATH"
RUN mkdir -p /work
WORKDIR /work
USER root
ENTRYPOINT []
CMD ["postgwas", "--help"]
