FROM ubuntu:22.04
MAINTAINER Melissa Meredith, melissa@datatecnica.com
RUN apt-get update && \
	DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata && \
	apt-get install -y cmake git make gcc g++ autoconf bzip2 wget unzip tabix libz-dev && \
	apt-get install -y samtools && \
	apt-get clean

# get hifiasm
ARG HIFIASM_VER=0.25.0
WORKDIR /opt
RUN wget https://github.com/chhylp123/hifiasm/archive/refs/tags/$HIFIASM_VER.zip && \
	unzip $HIFIASM_VER.zip && \
	rm $HIFIASM_VER.zip && \
	cd hifiasm-$HIFIASM_VER && \
	make 




WORKDIR /opt

