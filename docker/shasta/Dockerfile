FROM ubuntu:20.04
MAINTAINER Mikhail Kolmogorov, mkolmogo@ucsc.edu
RUN apt-get update && \
	DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata && \
	apt-get install -y cmake git make gcc g++ autoconf bzip2 wget tabix libz-dev && \
	apt-get install -y samtools && \
	apt-get clean

ARG SHASTA_VER=0.11.1
RUN wget https://github.com/paoloshasta/shasta/releases/download/$SHASTA_VER/shasta-Linux-$SHASTA_VER && \
	chmod +x shasta-Linux-$SHASTA_VER && \
	cp shasta-Linux-$SHASTA_VER /usr/bin/shasta

WORKDIR /opt
COPY ./shasta_config /opt/shasta_config
