FROM ubuntu

RUN apt -y update && apt install -y procps git python3 python3-pip

RUN mkdir -p /app/
WORKDIR /app/
RUN git clone https://gitlab.com/Sunmish/diffusefilter.git
WORKDIR diffusefilter
RUN pip install wheel
RUN python3 setup.py install 
WORKDIR scripts