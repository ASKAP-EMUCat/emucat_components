FROM ubuntu:focal

ARG DEBIAN_FRONTEND=noninteractive

RUN apt update
RUN apt -y upgrade
RUN apt install -y python3-setuptools
RUN apt install -y python3
RUN apt install -y wget
RUN apt install -y python3-pip
RUN pip3 install pyvo retry jinja2 asyncio asyncpg aiofiles

RUN mkdir /input
RUN mkdir /output
RUN mkdir /scripts
COPY ./scripts/casda.py /scripts
COPY ./scripts/catalog.py /scripts
WORKDIR /scripts
RUN chmod +x casda.py
RUN chmod +x catalog.py