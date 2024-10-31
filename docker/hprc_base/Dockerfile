FROM ubuntu:24.04
LABEL maintainer="wen-wei.liao@yale.edu"

USER root
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get install -y \
        git wget autoconf build-essential zlib1g-dev libbz2-dev libcurl4-gnutls-dev \
	liblzma-dev libncurses5-dev libncursesw5-dev libssl-dev

# miniconda3 (python v3.9)
ENV PATH="/opt/conda/bin:$PATH"
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py39_24.9.2-0-Linux-x86_64.sh -O ~/miniconda.sh && \
    bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    echo '. /opt/conda/etc/profile.d/conda.sh' > /etc/profile.d/conda.sh

# cmake v3.30.5
WORKDIR /opt/cmake_install
RUN mkdir /opt/cmake && \
    wget https://github.com/Kitware/CMake/releases/download/v3.30.5/cmake-3.30.5-linux-x86_64.sh && \
    bash cmake-3.30.5-linux-x86_64.sh --prefix=/opt/cmake --skip-license && \
    ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake && \
    rm -r /opt/cmake_install

# htslib v1.21
WORKDIR /opt/htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2 && \
    tar jxf htslib-1.21.tar.bz2 && \
    rm htslib-1.21.tar.bz2 && \
    cd htslib-1.21 && \
    ./configure && \
    make && \
    make install

# samtools v1.21
WORKDIR /opt/samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 && \
    tar jxf samtools-1.21.tar.bz2 && \
    rm samtools-1.21.tar.bz2 && \
    cd samtools-1.21 && \
    autoheader && \
    autoconf -Wno-syntax && \
    ./configure && \
    make && \
    make install

# bcftools v1.21
WORKDIR /opt/bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2 && \
    tar jxf bcftools-1.21.tar.bz2 && \
    rm bcftools-1.21.tar.bz2 && \
    cd bcftools-1.21/ && \
    autoheader && \
    autoconf -Wno-syntax && \
    ./configure && \
    make && \
    make install

WORKDIR /data
