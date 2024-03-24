FROM ubuntu:20.04
MAINTAINER Mikhail Kolmogorov, mkolmogo@ucsc.edu
RUN apt-get update && \
	DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata && \
	apt-get install -y cmake git make gcc g++ autoconf bzip2 wget tabix libz-dev libncurses5-dev libbz2-dev liblzma-dev && \
	DEBIAN_FRONTEND="noninteractive" apt-get install -y python3.10 && \
	apt-get clean

#RUN pip install sniffles
RUN update-alternatives --install /usr/bin/python3 python /usr/bin/python3.10 2
RUN update-alternatives --install /usr/bin/python3 python /usr/bin/python3.8 1

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


WORKDIR /opt/sniffles
RUN git clone https://github.com/fritzsedlazeck/Sniffles.git && \
    mv /opt/sniffles/Sniffles/src/sniffles/sniffles /opt/sniffles/Sniffles/src/sniffles2 && \
    touch /opt/sniffles/Sniffles/src/utils/__init__.py

ENV PATH="/opt/sniffles/Sniffles/src:${PATH}"

RUN mkdir -p /data
WORKDIR /data