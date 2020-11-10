FROM nvidia/cuda:10.2-devel

RUN apt update
RUN apt install -y wget git cmake python3 python3-pip build-essential libssl-dev libffi-dev python3-dev
RUN pip3 install --upgrade pip
WORKDIR /
RUN git clone https://github.com/HITS-AIN/PINK.git
WORKDIR PINK
RUN cmake .
RUN make install
WORKDIR /
RUN git clone https://github.com/tjgalvin/pyink.git
WORKDIR pyink
RUN pip3 install astropy extension-helpers numpy jinja2 sklearn scipy
RUN python3 setup.py install
RUN pip3 install pandas
WORKDIR /
RUN git clone https://github.com/avantyghem/ComplexXID.git