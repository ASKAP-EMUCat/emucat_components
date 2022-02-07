FROM ubuntu:focal

ARG DEBIAN_FRONTEND=noninteractive

RUN apt update
RUN apt -y upgrade
RUN apt install -y python3-setuptools
RUN apt install -y python3
RUN apt install -y wget
RUN apt install -y python3-pip
RUN pip3 install setuptools wheel numpy astropy scipy matplotlib pandas pyregion

RUN mkdir /input
RUN mkdir /output
RUN mkdir /scripts
COPY ./likelihood_ratio/likelihood_ratio_matching.py /scripts
COPY ./likelihood_ratio/lr_wrapper_emucat.py /scripts
WORKDIR /scripts
RUN chmod +x lr_wrapper_emucat.py