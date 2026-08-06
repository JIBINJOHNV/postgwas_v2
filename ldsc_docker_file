# LDSC Dockerfile (Python 2.7, amd64)
FROM --platform=linux/amd64 continuumio/miniconda3:4.7.12

WORKDIR /opt/ldsc

# 1) Get LDSC source
RUN conda install -y -n base -c conda-forge git && \
    conda clean -afy && \
    git clone https://github.com/bulik/ldsc.git .

# 2) Configure conda channels + create Python 2.7 env for LDSC
RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda config --set channel_priority strict && \
    conda create -y -n ldsc \
        python=2.7 \
        numpy=1.16 \
        scipy=1.2.1 \
        pandas=0.24.2 \
        bitarray \
        nose \
        pybedtools && \
    conda clean -afy

# 3) Put ldsc env + repo on PATH so ldsc.py / munge_sumstats.py are runnable
ENV PATH="/opt/conda/envs/ldsc/bin:/opt/ldsc:${PATH}"

# 4) Optional: smoke-test at build time
RUN /opt/conda/envs/ldsc/bin/python ldsc.py -h && \
    /opt/conda/envs/ldsc/bin/python munge_sumstats.py -h

# 5) Entry-point script that activates env then runs your command
COPY ldsc-entrypoint.sh /usr/local/bin/ldsc-entrypoint.sh
RUN chmod +x /usr/local/bin/ldsc-entrypoint.sh

ENTRYPOINT ["/usr/local/bin/ldsc-entrypoint.sh"]
CMD []
