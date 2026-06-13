FROM rocker/r-ver:4.6.0

ENV DEBIAN_FRONTEND=noninteractive
ENV RUSTUP_HOME=/usr/local/rustup
ENV CARGO_HOME=/usr/local/cargo
ENV PATH=/usr/local/cargo/bin:${PATH}

RUN apt-get update && apt-get install -y --no-install-recommends \
    bash \
    ca-certificates \
    curl \
    g++ \
    gfortran \
    git \
    libbz2-dev \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libfribidi-dev \
    libgit2-dev \
    libharfbuzz-dev \
    libjpeg-dev \
    liblzma-dev \
    libmagick++-dev \
    libpng-dev \
    libssl-dev \
    libtiff-dev \
    libxml2-dev \
    make \
    pandoc \
    pkg-config \
    qpdf \
    zlib1g-dev \
  && rm -rf /var/lib/apt/lists/*

RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \
    | sh -s -- -y --profile minimal --default-toolchain stable \
  && rustc --version \
  && cargo --version

WORKDIR /workspace

COPY DESCRIPTION /workspace/DESCRIPTION
COPY scripts/flounder-dependencies.R /workspace/scripts/flounder-dependencies.R
COPY scripts/bootstrap-r-dependencies.R /workspace/scripts/bootstrap-r-dependencies.R
COPY scripts/audit-r-dependencies.R /workspace/scripts/audit-r-dependencies.R

RUN Rscript scripts/bootstrap-r-dependencies.R --install \
  && Rscript scripts/audit-r-dependencies.R --output=/tmp/flounder-r-dependencies.tsv

COPY . /workspace

CMD ["bash"]
