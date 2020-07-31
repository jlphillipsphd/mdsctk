FROM debian:buster

LABEL maintainer="Joshua L. Phillips <https://www.cs.mtsu.edu/~jphillips/>"

USER root

# Dependencies in APT repos
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -yq \
    build-essential \
    cmake \
    git \
    libarpack++2-dev \
    libboost-program-options-dev \
    libdb++-dev \
    libgromacs-dev \
    pkg-config \
    python \
    r-base \
    && apt-get clean

# Dependencies -not- in APT repos
WORKDIR /root
RUN R -e 'install.packages("argparse")'
RUN R -e 'install.packages("fields")'

# Code build
WORKDIR /opt
RUN git clone https://github.com/jlphillipsphd/mdsctk.git
WORKDIR /opt/mdsctk
RUN git checkout development
WORKDIR /opt/mdsctk
RUN cmake . && make

# Deployment
WORKDIR /opt/mdsctk/examples
ENV MDSCTK_HOME="/opt/mdsctk"
ENV PATH="$PATH:/opt/mdsctk"
