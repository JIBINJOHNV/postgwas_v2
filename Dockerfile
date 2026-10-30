# Ubuntu and native installations use the same installer and version manifest.
# Build locally with --platform linux/amd64; review MAGMA licensing before
# redistributing an image containing its binary.
FROM ubuntu:24.04@sha256:fe81560990fae86dc4743202b2cb4d3df98d0bcf949173cecef15ce1127ef852
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates curl git bzip2 unzip make gcc g++ file procps && \
    rm -rf /var/lib/apt/lists/*

ARG MINIFORGE_VERSION=26.9.4-1
ARG MINIFORGE_SHA256=571b410272cc434f28bc6c5c49ca8ffd4a4d1bcba9f57f1b6baccdb46446cfcc
RUN curl --fail --location --silent --show-error \
    "https://github.com/conda-forge/miniforge/releases/download/${MINIFORGE_VERSION}/Miniforge3-${MINIFORGE_VERSION}-Linux-x86_64.sh" \
    -o /tmp/miniforge.sh && \
    echo "${MINIFORGE_SHA256}  /tmp/miniforge.sh" | sha256sum --check - && \
    bash /tmp/miniforge.sh -b -p /opt/conda && rm /tmp/miniforge.sh
ENV PATH="/opt/conda/bin:${PATH}"

WORKDIR /opt/postgwas
# Copy only install inputs; study data and local notes do not enter the image.
COPY pyproject.toml README.md environment.yml ./
COPY src ./src
COPY tools ./tools
ARG POSTGWAS_INSTALL_JOBS=2
RUN bash tools/setup/install_postgwas.sh --all-tools --name postgwas \
    --jobs "${POSTGWAS_INSTALL_JOBS}"

# Conda activation supplies the same isolated runtime settings as native use.
# A writable HOME also permits an explicitly selected non-root Docker user.
ENV HOME=/tmp
WORKDIR /work
ENTRYPOINT ["/opt/conda/bin/conda", "run", "--no-capture-output", "-n", "postgwas"]
CMD ["postgwas", "--help"]
