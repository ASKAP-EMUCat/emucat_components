FROM ubuntu:focal

ARG DEBIAN_FRONTEND=noninteractive

RUN apt update
RUN apt -y upgrade
RUN apt install -y python3-setuptools
RUN apt install -y python3
RUN apt install -y wget
RUN apt install -y python3-pip
RUN apt install -y libcurl4-openssl-dev libssl-dev
RUN pip3 install wheel
RUN pip3 install astropy pyvo retry jinja2 asyncio asyncpg aiofiles
RUN pip3 install --ignore-installed --no-cache-dir noaodatalab

RUN mkdir /input
RUN mkdir /output
RUN mkdir /scripts
COPY ./scripts/casda.py /scripts
COPY ./scripts/catalog.py /scripts
COPY ./scripts/properties.py /scripts
COPY ./scripts/noao.py /scripts
WORKDIR /scripts
RUN chmod +x casda.py
RUN chmod +x catalog.py
RUN chmod +x noao.py