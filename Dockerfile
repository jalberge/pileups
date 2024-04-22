FROM gcr.io/broad-getzlab-workflows/base_image:v0.0.6

ARG VERSION=1.20
WORKDIR /usr/bin

RUN apt-get update && apt-get install -y --no-install-recommends \
        autoconf \
        automake \
        gcc \
        perl \
        bzip2 \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-gnutls-dev \
        libssl-dev \
        libncurses5-dev && \
        rm -rf /var/lib/apt/lists/*

RUN  curl -L -O https://github.com/samtools/htslib/releases/download/${VERSION}/htslib-${VERSION}.tar.bz2 \
        && tar -vxjf htslib-${VERSION}.tar.bz2 \
        && cd htslib-${VERSION} \
        && make

RUN curl -L -O https://github.com/samtools/samtools/releases/download/${VERSION}/samtools-${VERSION}.tar.bz2 \
        && tar -vxjf samtools-${VERSION}.tar.bz2 \
        && cd samtools-${VERSION} \
        && ./configure --enable-libcurl \
        && make all all-htslib \
        && make install install-htslib

RUN curl -L -O https://github.com/samtools/bcftools/releases/download/${VERSION}/bcftools-${VERSION}.tar.bz2 \
        && tar -vxjf bcftools-${VERSION}.tar.bz2 \
        && cd bcftools-${VERSION} \
        && ./configure --enable-libcurl \
        && make all all-htslib \
        && make install install-htslib

ENV PATH="/usr/bin/bcftools-${VERSION}:${PATH}"
ENV PATH="/usr/bin/samtools-${VERSION}:${PATH}"
ENV PATH="/usr/bin/htslib-${VERSION}:${PATH}"
