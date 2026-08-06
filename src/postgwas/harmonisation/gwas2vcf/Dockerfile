# Use a small base image
FROM python:3.8-slim

LABEL maintainer="Jibin John: johnjibinv@gmail.com" \
      description="GWAS sums stat harmonisaiton and for converting to vcf format and liftover;  with CrossMap, vgraph, bcftools, polars, pandas, pyarrow"

# Install only what's necessary
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        bcftools \
        tabix \
        git \
        build-essential \
        libbz2-dev \
        liblzma-dev \
        libz-dev \
        libcurl4 \
        ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app
COPY . /app

# Install Python dependencies in one layer
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir pandas polars pyarrow scipy && \
    pip install --no-cache-dir -r requirements.txt && \
    pip install --no-cache-dir git+https://github.com/bioinformed/vgraph@v1.4.0#egg=vgraph && \
    pip install --no-cache-dir git+https://github.com/liguowang/CrossMap.git && \
    # remove compiler tools to reduce final image size
    apt-get purge -y build-essential libbz2-dev liblzma-dev libz-dev && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/* /root/.cache/pip

# Optional Flask app port
EXPOSE 5000

# Default command
CMD ["python", "main.py"]

