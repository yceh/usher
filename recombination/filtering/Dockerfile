FROM ubuntu:20.04
ENV APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE=DontWarn
ENV DEBIAN_FRONTEND=noninteractive
USER root

RUN apt-get update && apt-get install -yq --no-install-recommends \
    build-essential git wget vim curl rsync python3 python3-pip cmake ninja-build jq \
    bzip2 gnupg2 squashfs-tools openmpi-bin \
    libboost-all-dev \
    libprotoc-dev libprotoc-dev protobuf-compiler \
    libtbb-dev \
    mpich libmpich-dev automake libtool autoconf make nasm \
    ca-certificates \
    apt-transport-https gnupg \
    lsb-core \
    sudo 

# gcloud
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
   curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - && \
   apt-get update -y && \
   apt-get install -y google-cloud-sdk

# gcsfuse
#ENV GCSFUSE_REPO=gcsfuse-`lsb_release -c -s`
ENV GCSFUSE_REPO=gcsfuse-focal
RUN echo "deb http://packages.cloud.google.com/apt $GCSFUSE_REPO main" | sudo tee /etc/apt/sources.list.d/gcsfuse.list && \
  curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | sudo apt-key add - && \
  apt-get update -y && \
  apt-get install -yq gcsfuse

# mafft build
RUN git clone https://github.com/GSLBiotech/mafft && \
    cd mafft/core && \
    make -j$(grep -c ^processor /proc/cpuinfo 2>/dev/null || sysctl -n hw.ncpu) && \
    make install

# install conda
RUN curl -Ol https://repo.anaconda.com/miniconda/Miniconda3-py38_4.10.3-Linux-x86_64.sh
RUN bash Miniconda3-py38_4.10.3-Linux-x86_64.sh -b

ENV PATH="/root/miniconda3/bin:${PATH}"

RUN conda install mamba -n base -c conda-forge

RUN mamba install -y -c conda-forge -c bioconda snakemake-minimal numpy pyyaml

# install faToVcf
RUN rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/faToVcf /usr/bin
RUN chmod +x /usr/bin/faToVcf

# install faSomeRecords
RUN rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/faSomeRecords /usr/bin
RUN chmod +x /usr/bin/faSomeRecords

# usher stuff
WORKDIR /HOME/kentsource
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faSomeRecords
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faSize

RUN chmod 775 *
WORKDIR /HOME

#RUN git clone https://github.com/yatisht/usher.git
#WORKDIR usher
RUN git clone https://github.com/mrkylesmith/usher.git
WORKDIR usher

RUN mkdir -p build && \
   cd build && \
   wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz && \
   tar -xvzf 2019_U9.tar.gz 

RUN mkdir -p build && \
   cd build && \
   wget https://github.com/intel/isa-l/archive/refs/tags/v2.30.0.tar.gz && \
   tar -xvf v2.30.0.tar.gz && \
   cd isa-l-2.30.0 && \
   ./autogen.sh && \
   ./configure && \
   make -j$(nproc) && \
   make install

RUN cd build && \
   cmake -GNinja -DTBB_DIR=${PWD}/oneTBB-2019_U9 -DCMAKE_PREFIX_PATH=${PWD}/oneTBB-2019_U9/cmake .. 

RUN cd build && ninja

# create top level data directory for mount, local file save option
RUN mkdir /data

#COPY process.py /data
#COPY run.py /data

# set the path
ENV PATH="/HOME/usher/build:/HOME/kentsource:${PATH}"
