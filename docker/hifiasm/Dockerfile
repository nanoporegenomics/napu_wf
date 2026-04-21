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

ARG RASUSA_VER=3.0.0
RUN wget https://github.com/mbhall88/rasusa/releases/download/3.0.0/rasusa-$RASUSA_VER-x86_64-unknown-linux-gnu.tar.gz && \
	tar -xzvf rasusa-$RASUSA_VER-x86_64-unknown-linux-gnu.tar.gz && \
	rm rasusa-$RASUSA_VER-x86_64-unknown-linux-gnu.tar.gz

ARG RASUSA_VER=3.0.0
RUN wget https://github.com/mbhall88/rasusa/releases/download/$RASUSA_VER/rasusa-$RASUSA_VER-x86_64-unknown-linux-gnu.tar.gz && \
    tar -xzvf rasusa-$RASUSA_VER-x86_64-unknown-linux-gnu.tar.gz && \
    mv rasusa-$RASUSA_VER-x86_64-unknown-linux-gnu/rasusa /usr/local/bin/rasusa && \
    chmod +x /usr/local/bin/rasusa && \
    rm -rf rasusa-$RASUSA_VER-x86_64-unknown-linux-gnu*



ENV PATH="/opt/hifiasm-$HIFIASM_VER:$PATH"

WORKDIR /data

