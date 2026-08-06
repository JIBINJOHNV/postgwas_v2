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
# 2. VEP & LDSC ENV (Combined Cleanups)
# =====================================================================
RUN micromamba create -y -n vep -c conda-forge -c bioconda ensembl-vep=113 && \
    micromamba clean --all --yes

RUN micromamba create -y -n ldsc -c conda-forge -c bioconda \
    python=2.7 numpy=1.16 scipy=1.2.1 pandas=0.24.2 bitarray nose pybedtools && \
    micromamba clean --all --yes

WORKDIR /opt
RUN git clone https://github.com/bulik/ldsc.git && \
    rm -rf ldsc/.git

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
RUN wget https://github.com/samtools/bcftools/releases/download/1.23/bcftools-1.23.tar.bz2 && \
    wget https://github.com/samtools/htslib/releases/download/1.23/htslib-1.23.tar.bz2 && \
    tar -xvjf bcftools-1.23.tar.bz2 && \
    tar -xvjf htslib-1.23.tar.bz2 && \
    cd htslib-1.23 && \
    ./configure --enable-libcurl --prefix=/usr/local && \
    make -j && make install && \
    cd ../bcftools-1.23/plugins && \
    wget https://raw.githubusercontent.com/freeseek/score/master/score.c && \
    wget https://raw.githubusercontent.com/freeseek/score/master/score.h && \
    wget https://raw.githubusercontent.com/freeseek/score/master/munge.c && \
    wget https://raw.githubusercontent.com/freeseek/score/master/liftover.c && \
    wget https://raw.githubusercontent.com/freeseek/score/master/metal.c && \
    wget https://raw.githubusercontent.com/freeseek/score/master/blup.c && \
    wget https://raw.githubusercontent.com/freeseek/score/master/pgs.c && \
    wget https://raw.githubusercontent.com/freeseek/score/master/pgs.mk && \
    cd .. && \
    sed -i '2254s/^/\/\//' plugins/pgs.c && \
    sed -i '2255s/^/\/\//' plugins/pgs.c && \
    ./configure --prefix=/usr/local --with-htslib=/opt/tools/htslib-1.23 CPPFLAGS="-I/usr/include/suitesparse" CFLAGS="-I/usr/include/suitesparse" && \
    make -j && make install && \
    cd /opt/tools && \
    rm -rf bcftools-1.23 htslib-1.23 *.tar.bz2

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

COPY magma/magma /usr/local/bin/magma
RUN chmod +x /usr/local/bin/magma

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
