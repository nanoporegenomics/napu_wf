FROM ubuntu:22.04

MAINTAINER Mikhail Kolmogorov, mkolmogo@ucsc.edu
RUN apt-get update && \
	DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata && \
	apt-get install -y cmake git make automake gcc g++ autoconf bzip2 \
                       wget tabix libz-dev libncurses5-dev libbz2-dev \
                       perl zlib1g-dev libcurl4-gnutls-dev liblzma-dev \
                       libssl-dev libperl-dev libgsl0-dev && \
	apt-get -y install bedtools && \
    apt-get -y install python3-pip && \
    apt-get clean


#RUN pip install sniffles
RUN python3 --version
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install pysam scipy edlib matplotlib biopython

WORKDIR /opt
COPY ./trf_annotations /opt/trf_annotations

### samtools
# 1.9
WORKDIR /opt/samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar xvf samtools-1.9.tar.bz2 && \
	rm -r /opt/samtools/samtools-1.9.tar.bz2 && \
	cd samtools-1.9/ && \
	autoheader && \
	autoconf -Wno-header && \
	./configure && \
	make && \
	cp samtools /usr/bin/samtools

ARG bcftoolsVer="1.19"
WORKDIR /opt/bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/${bcftoolsVer}/bcftools-${bcftoolsVer}.tar.bz2 && \
   tar -vxjf bcftools-${bcftoolsVer}.tar.bz2 && \
   rm bcftools-${bcftoolsVer}.tar.bz2 && \
   cd bcftools-${bcftoolsVer} && \
   make && \
   make install

WORKDIR /opt/sniffles
RUN git clone https://github.com/fritzsedlazeck/Sniffles.git && \
    mv /opt/sniffles/Sniffles/src/sniffles/sniffles /opt/sniffles/Sniffles/src/sniffles2 && \
    touch /opt/sniffles/Sniffles/src/utils/__init__.py

ENV PATH="/opt/sniffles/Sniffles/src:${PATH}"
ENV PATH="/opt/bcftools/bcftools2:${PATH}"

RUN mkdir -p /data
WORKDIR /data