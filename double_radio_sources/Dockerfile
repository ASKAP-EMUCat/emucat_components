FROM ubuntu:focal

ARG DEBIAN_FRONTEND=noninteractive

RUN apt update
RUN apt -y upgrade
RUN apt install -y python3-setuptools
RUN apt install -y python3
RUN apt install -y wget
RUN apt install -y python3-pip
RUN pip3 install wheel numpy astropy numpy argparse pandas

RUN mkdir /input
RUN mkdir /output
RUN mkdir /scripts
COPY ./double_radio_sources/emu_doubles.py /scripts
WORKDIR /scripts
RUN chmod +x emu_doubles.py