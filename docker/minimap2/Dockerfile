FROM ubuntu:20.04
MAINTAINER Mikhail Kolmogorov, mkolmogo@ucsc.edu
RUN apt-get update && \
	DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata && \
	apt-get install -y cmake git make gcc g++ autoconf bzip2 wget tabix libz-dev libncurses5-dev libbz2-dev liblzma-dev && \
	apt-get clean

ARG MM_VER=2.23
RUN wget https://github.com/lh3/minimap2/releases/download/v$MM_VER/minimap2-$MM_VER.tar.bz2 && \
	tar xvf minimap2-$MM_VER.tar.bz2 && \
	rm minimap2-$MM_VER.tar.bz2 && \
	cd minimap2-$MM_VER && \
	make && \
	cp minimap2 /usr/bin/

### samtools
# 1.15
WORKDIR /opt/samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.15/samtools-1.15.tar.bz2 && \
    tar xvf samtools-1.15.tar.bz2 && \
	rm -r /opt/samtools/samtools-1.15.tar.bz2 && \
	cd samtools-1.15/ && \
	autoheader && \
	autoconf -Wno-header && \
	./configure && \
	make && \
	cp samtools /usr/bin/samtools
