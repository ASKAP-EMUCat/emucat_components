FROM ubuntu:focal

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Australia
ARG COMPILE_JOBS=2
ARG YANDASOFT_BRANCH=develop

RUN apt-get update

# Toolchains, etc
RUN apt-get install -y \
    git \
    cmake \
    flex \
    bison \
    gfortran \
    git \
    g++ \
    make \
    patch \
    subversion \
    wget

# Casacore/yandasoft dependencies
RUN apt-get install -y \
    casacore-dev \
    libboost-dev \
    libboost-all-dev \
    libcfitsio-dev \
    libcppunit-dev \
    libfftw3-dev \
    libgsl-dev \
    liblog4cxx-dev \
    libopenblas-dev \
    libopenmpi-dev \
    libxerces-c-dev \
    wcslib-dev \
    libopenmpi3

#############################################################
# Build casacore, then all yandasoft dependencies, then yandasoft itself
WORKDIR /home
RUN git clone https://github.com/casacore/casarest.git
WORKDIR /home/casarest
RUN git checkout v1.7.0
RUN mkdir build
WORKDIR build
RUN cmake ..
RUN cmake --build . -j $COMPILE_JOBS
RUN cmake --build . --target install -j $COMPILE_JOBS

WORKDIR /home
RUN git clone https://github.com/rtobar/all_yandasoft.git
WORKDIR /home/all_yandasoft
RUN ./git-do clone -b $YANDASOFT_BRANCH
RUN mkdir build
WORKDIR build
RUN cmake -DENABLE_OPENMP=YES ..
RUN cmake --build . -j $COMPILE_JOBS
RUN cmake --build . --target install -j $COMPILE_JOBS
