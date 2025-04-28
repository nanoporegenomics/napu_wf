FROM ubuntu:22.04
MAINTAINER Melissa Meredith, melissa@datatecnica.com

# dependencies
RUN apt-get update && \
    DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata && \
    apt-get install -y cmake git make gcc g++ autoconf bzip2 wget unzip tabix libz-dev \
        curl python3 python3-pip ca-certificates apt-transport-https lsb-release gnupg && \
    apt-get install -y samtools && \
    apt-get clean


# install the Google Cloud SDK (includes gsutil)
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" \
    | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg \
    | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add - && \
    apt-get update && apt-get install -y google-cloud-sdk


