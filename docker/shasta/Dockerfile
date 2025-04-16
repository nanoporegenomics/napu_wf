FROM ubuntu:20.04
MAINTAINER Mikhail Kolmogorov, mkolmogo@ucsc.edu
RUN apt-get update && \
	DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata && \
	apt-get install -y cmake git make gcc g++ autoconf bzip2 wget tabix libz-dev && \
	apt-get install -y samtools && \
	apt-get clean

ARG SHASTA_VER=0.14.0
RUN wget https://github.com/paoloshasta/shasta/releases/download/$SHASTA_VER/shasta-Linux-$SHASTA_VER && \
	chmod +x shasta-Linux-$SHASTA_VER && \
	cp shasta-Linux-$SHASTA_VER /usr/bin/shasta

# chopper for Qscore cutoffs
ARG CHOPPER_VER=0.9.2
RUN wget https://github.com/wdecoster/chopper/releases/download/v$CHOPPER_VER/chopper-linux.zip && \
	unzip chopper-linux.zip && \
	chmod +x chopper
	cp chopper /usr/bin/chopper


WORKDIR /opt
COPY ./shasta_config /opt/shasta_config
