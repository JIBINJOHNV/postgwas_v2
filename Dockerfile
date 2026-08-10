# =====================================================================
# MiXeR build stage
# =====================================================================
# This temporary stage compiles MiXeR. Boost, source code, Git history,
# tests and build files are not copied into the final PostGWAS image.
FROM mambaorg/micromamba:1.4.2 AS mixer-builder

USER root

ARG GSA_MIXER_TAG=v2.2.1
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
# required by MiXeR's FastDifferentialCoding implementation.
RUN git clone \
    --depth 1 \
    --branch "${GSA_MIXER_TAG}" \
    https://github.com/precimed/gsa-mixer.git \
    /tmp/gsa-mixer && \
    sed -i 's/-march=native/-mssse3 -msse4.1/g' \
    /tmp/gsa-mixer/src/CMakeLists.txt && \
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
FROM mambaorg/micromamba:1.4.2

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
    libtiff5-dev \
    libjpeg-dev \
    libfontconfig1-dev \
    libgit2-dev \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libeigen3-dev \
    procps && \
    rm -rf /var/lib/apt/lists/*

# =====================================================================
# 1. MAIN ENV: PostGWAS
# =====================================================================
COPY environment.yml /tmp/environment.yml

# Create env and clean immediately.
RUN micromamba create -y -n postgwas -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Install Build Tools -> Install R/Deps -> Install minimal MiXeR Python
# dependencies -> CLEAN UP.
RUN micromamba install -y -n postgwas -c conda-forge \
    cmake \
    mkl \
    mkl-include \
    sysroot_linux-64 \
    gcc_linux-64 \
    gxx_linux-64 \
    r-base \
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
    six && \
    micromamba run -n postgwas Rscript -e "install.packages('remotes', repos='https://cloud.r-project.org'); remotes::install_github('jianyanglab/gsmr2', dependencies=TRUE)" && \
    micromamba clean --all --yes

ENV PATH="/opt/conda/envs/postgwas/bin:$PATH"

# =====================================================================
# 2. VEP
# =====================================================================
RUN micromamba create -y -n vep -c conda-forge -c bioconda ensembl-vep=113 && \
    micromamba clean --all --yes

# LDSC runs in the main PostGWAS environment. Pinning the upstream commit keeps
# image builds reproducible while using CBIIT's maintained Python 3 branch.
ARG LDSC_REPOSITORY=https://github.com/CBIIT/ldsc.git
ARG LDSC_COMMIT=6c673952cee74bd5c57aef1555a03b1c015399a0
RUN mkdir -p /opt/ldsc && \
    git -C /opt/ldsc init && \
    git -C /opt/ldsc remote add origin "${LDSC_REPOSITORY}" && \
    git -C /opt/ldsc fetch --depth 1 origin "${LDSC_COMMIT}" && \
    git -C /opt/ldsc checkout --detach FETCH_HEAD && \
    test "$(git -C /opt/ldsc rev-parse HEAD)" = "${LDSC_COMMIT}" && \
    rm -rf /opt/ldsc/.git && \
    micromamba run -n postgwas \
    python -m pip install --no-deps --no-cache-dir /opt/ldsc && \
    micromamba run -n postgwas ldsc.py --help >/dev/null && \
    micromamba run -n postgwas munge_sumstats.py --help >/dev/null && \
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
# 3. ENRICHER ENV (New Module)
# =====================================================================
# We combine the requested packages into one create command.
# Note: rpy2 will automatically pull a compatible r-base into this env.
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
    rich-argparse \
    scipy \
    "numpy<2.0" \
    r-webgestaltr \
    zip \
    gprofiler-official && \
    micromamba clean --all --yes

# =====================================================================
# Build HTSlib + BCFtools (And DELETE SOURCE afterwards)
# =====================================================================
WORKDIR /opt/tools
SHELL ["/bin/bash", "-c"]

# Download -> Compile -> Install -> DELETE SOURCE.
RUN wget https://github.com/samtools/bcftools/releases/download/1.23.1/bcftools-1.23.1.tar.bz2 && \
    wget https://github.com/samtools/htslib/releases/download/1.23.1/htslib-1.23.1.tar.bz2 && \
    tar -xvjf bcftools-1.23.1.tar.bz2 && \
    tar -xvjf htslib-1.23.1.tar.bz2 && \
    cd htslib-1.23.1 && \
    ./configure --enable-libcurl --prefix=/usr/local && \
    make -j && make install && \
    cd ../bcftools-1.23.1/plugins && \
    wget https://raw.githubusercontent.com/freeseek/score/909d23019e19aeadf3bf6fe1407fd6afc094592a/score.c && \
    wget https://raw.githubusercontent.com/freeseek/score/909d23019e19aeadf3bf6fe1407fd6afc094592a/score.h && \
    wget https://raw.githubusercontent.com/freeseek/score/909d23019e19aeadf3bf6fe1407fd6afc094592a/munge.c && \
    wget https://raw.githubusercontent.com/freeseek/score/909d23019e19aeadf3bf6fe1407fd6afc094592a/liftover.c && \
    wget https://raw.githubusercontent.com/freeseek/score/909d23019e19aeadf3bf6fe1407fd6afc094592a/metal.c && \
    wget https://raw.githubusercontent.com/freeseek/score/909d23019e19aeadf3bf6fe1407fd6afc094592a/blup.c && \
    wget https://raw.githubusercontent.com/freeseek/score/909d23019e19aeadf3bf6fe1407fd6afc094592a/pgs.c && \
    wget https://raw.githubusercontent.com/freeseek/score/909d23019e19aeadf3bf6fe1407fd6afc094592a/pgs.mk && \
    cd .. && \
    sed -i '2254s/^/\/\//' plugins/pgs.c && \
    sed -i '2255s/^/\/\//' plugins/pgs.c && \
    ./configure --prefix=/usr/local --with-htslib=/opt/tools/htslib-1.23.1 CPPFLAGS="-I/usr/include/suitesparse" CFLAGS="-I/usr/include/suitesparse" && \
    make -j && make install && \
    cd /opt/tools && \
    rm -rf bcftools-1.23.1 htslib-1.23.1 *.tar.bz2

ENV BCFTOOLS_PLUGINS="/usr/local/libexec/bcftools"

# =====================================================================
# Install LDStore, FINEMAP, GCTA, MAGMA
# =====================================================================
WORKDIR /tmp/install_tools

RUN wget http://www.christianbenner.com/ldstore_v2.0_x86_64.tgz && \
    tar -xzf ldstore_v2.0_x86_64.tgz && \
    mv ldstore_v2.0_x86_64/ldstore_v2.0_x86_64 /usr/local/bin/ldstore && \
    chmod +x /usr/local/bin/ldstore && \
    wget http://www.christianbenner.com/finemap_v1.4.2_x86_64.tgz && \
    tar -xzf finemap_v1.4.2_x86_64.tgz && \
    mv finemap_v1.4.2_x86_64/finemap_v1.4.2_x86_64 /usr/local/bin/finemap && \
    chmod +x /usr/local/bin/finemap && \
    wget --user-agent="Mozilla/5.0" https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.95.0-linux-kernel-3-x86_64.zip && \
    unzip gcta-1.95.0-linux-kernel-3-x86_64.zip && \
    mv gcta-1.95.0-linux-kernel-3-x86_64/gcta64 /usr/local/bin/gcta && \
    chmod +x /usr/local/bin/gcta && \
    wget https://yanglab.westlake.edu.cn/software/smr/download/smr-1.3.1-linux-x86_64.zip && \
    unzip smr-1.3.1-linux-x86_64.zip && \
    mv smr-1.3.1-linux-x86_64/smr /usr/local/bin/smr && \
    chmod +x /usr/local/bin/smr && \
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
# Install SnpEff
# =====================================================================
WORKDIR /opt

RUN wget https://snpeff-public.s3.amazonaws.com/versions/snpEff_latest_core.zip && \
    unzip snpEff_latest_core.zip && \
    rm snpEff_latest_core.zip && \
    rm -rf /opt/snpEff/examples

ENV SNPEFF_HOME=/opt/snpEff
ENV PATH="${SNPEFF_HOME}:$PATH"

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
# Install PostGWAS (Keep R-related dependencies)
# =====================================================================
WORKDIR /opt/postgwas
COPY . /opt/postgwas

# ONLY remove heavy build-only tools, NOT the sysroot or compilers
# that R and its libraries depend on.
RUN micromamba run -n postgwas pip install --upgrade pip && \
    micromamba run -n postgwas pip install --no-deps --no-cache-dir -e . && \
    micromamba run -n postgwas Rscript -e \
    "parse(file='/opt/postgwas/src/postgwas/modules/caldera/run_caldera.R')" \
    >/dev/null && \
    micromamba remove -n postgwas -y \
    cmake \
    mkl-include && \
    micromamba clean --all --yes

# =====================================================================
# Final Config
# =====================================================================
ENV HOME=/tmp
RUN chown -R mambauser:mambauser /opt/postgwas
ENV PATH="/opt/conda/envs/postgwas/bin:$PATH"
USER root
ENTRYPOINT []
