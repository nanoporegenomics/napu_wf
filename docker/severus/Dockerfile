FROM ubuntu:22.04
MAINTAINER Melissa Meredith, melissa@datatecnica.com

RUN apt-get update && \
	DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata && \
	apt-get install -y cmake git make gcc g++ autoconf bzip2 wget unzip tabix libz-dev && \
	apt-get install -y samtools ca-certificates && \
	apt-get clean

# install miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
	bash miniconda.sh -b -p /opt/conda && \
	rm miniconda.sh

# add conda to PATH
ENV PATH=/opt/conda/bin:$PATH

# use bash for RUN commands from here
SHELL ["/bin/bash","-c"]

# accept channel ToS
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
	conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r 
	
# clone the Severus git
RUN git clone https://github.com/KolmogorovLab/Severus /opt/Severus


# create create conda environment
RUN conda env create --name severus_env --file /opt/Severus/environment.yml 

ENV PATH=/opt/conda/envs/severus_env/bin:$PATH

RUN echo "source /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc
