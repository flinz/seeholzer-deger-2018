# Pull base image.
FROM ubuntu:16.04

MAINTAINER "Alex Seeholzer" <seeholzer@gmail.com>

ENV TERM=xterm \
    TZ=Europe/Berlin \
    DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
        software-properties-common \
        tzdata \
        locales \
        ssh \
        curl \
        wget \
        git \
        unzip \
        build-essential \
        autoconf \
        libltdl-dev \
        libreadline6-dev \
        libncurses5-dev \
        libgsl-dev \
        openmpi-bin \
        libopenmpi-dev \
        libhdf5-dev \
        pep8 \
        python-minimal \
        python-dev \
        python-setuptools \
        python-path \
        python-numpy \
        python-matplotlib \
        python-scipy \
        python-pip \
        cython=0.23.4-0ubuntu5 && \
        apt-get autoremove -y && \
        rm -rf /var/lib/apt/lists/* && \
        locale-gen en_US.UTF-8

ENV LANG='en_US.UTF-8' LANGUAGE='en_US:en' LC_ALL='en_US.UTF-8'

RUN curl https://bootstrap.pypa.io/get-pip.py -o /tmp/get-pip.py; python /tmp/get-pip.py;

COPY Makefile /seeholzer-deger-2018/Makefile
COPY nest /seeholzer-deger-2018/nest
WORKDIR /seeholzer-deger-2018/
RUN make nest
COPY . /seeholzer-deger-2018
RUN HDF5_DIR=/usr/lib/x86_64-linux-gnu/hdf5/serial pip install -r requirements.txt

ENV PYTHONPATH="/seeholzer-deger-2018/nest/nest_build/lib/python2.7/site-packages:/seeholzer-deger-2018/nest/nest_build/lib64/python2.7/site-packages:${PYTHONPATH}"
ENV PATH="/seeholzer-deger-2018/nest/nest_build/bin:${PATH}"
ENV DYLD_LIBRARY_PATH="/seeholzer-deger-2018/nest/nest_build/lib/nest:${DYLD_LIBRARY_PATH}"
ENV LD_LIBRARY_PATH="/seeholzer-deger-2018/nest/nest_build/lib/nest:${LD_LIBRARY_PATH}"