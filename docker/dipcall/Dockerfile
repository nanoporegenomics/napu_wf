FROM ubuntu:20.04
MAINTAINER Mikhail Kolmogorov, mkolmogo@ucsc.edu
RUN apt-get update && \
	DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata && \
	apt-get install -y cmake git make gcc g++ autoconf bzip2 wget tabix libz-dev libncurses5-dev libbz2-dev liblzma-dev && \
	apt-get install -y python3-pip && \
	apt-get clean

RUN pip install sniffles

WORKDIR /opt
RUN wget https://github.com/lh3/dipcall/releases/download/v0.3/dipcall-0.3_x64-linux.tar.bz2 && tar -jxf dipcall-0.3_x64-linux.tar.bz2
