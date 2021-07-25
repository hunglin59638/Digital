FROM ubuntu

LABEL base.image="ubuntu"
LABEL dockerfile.version="1"
LABEL software="Digital"
LABEL software.version="0.0.0"
LABEL description="A pipeline for annotation of proteins sequence. It could annoate with eggnog-mapper and then alingn NR database with diamond."
LABEL website="https://github.com/hunglin59638/digitama"
# LABEL license="https://github.com/hunglin59638/digitama/blob/master/LICENSE"
LABEL maintainer="Hung-Lin Chen"
LABEL maintainer.email="hunglin59638@gmail.com"

RUN apt-get update && apt-get install -y \
 python \
 python3 \
 python3-pip \
 git \
 wget

RUN mkdir -p /opt/digital
WORKDIR /opt/digital
ADD . .
# COPY modules/* modules/
RUN python3 setup.py
WORKDIR lib
RUN wget -qO- http://github.com/bbuchfink/diamond/releases/download/v2.0.8/diamond-linux64.tar.gz | tar xzf -
RUN mkdir /opt/digital/db
RUN  mkdir /data
WORKDIR /data

ENV PATH="/opt/digital:/opt/digital/lib/eggnog-mapper:/opt/digital/lib:${PATH}"