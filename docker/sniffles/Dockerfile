FROM ubuntu:22.04
MAINTAINER Melissa Meredith, melissa@datatecnica.com
RUN apt-get update && \
	DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata && \
	apt-get install -y cmake git make gcc g++ autoconf bzip2 wget tabix libz-dev libncurses5-dev libbz2-dev liblzma-dev && \
	apt-get install -y python3-pip && \
	apt-get install -y samtools && \
	apt-get clean

RUN pip install sniffles

WORKDIR /opt
COPY ./trf_annotations /opt/trf_annotations

