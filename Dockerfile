FROM nvidia/cuda:11.7.1-devel-ubuntu20.04

RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    gfortran \
    flex \
    bison \
    tcsh \
    csh \
    xorg-dev \
    libxmu-dev \
    libxi-dev \
    wget \
    git \
    python3 python3-pip \
    libnetcdf-dev \
    libopenmpi-dev openmpi-bin && \
    rm -rf /var/lib/apt/lists/*

#AmberTools22
WORKDIR /opt
RUN wget https://ambermd.org/downloads/AmberTools22.tar.bz2 && \
    tar -xjf AmberTools22.tar.bz2 && \
    rm AmberTools22.tar.bz2

WORKDIR /opt/amber22
RUN ./configure gnu && \
    make install -j$(nproc)

ENV AMBERHOME=/opt/amber22
ENV PATH="$AMBERHOME/bin:$PATH"
ENV LD_LIBRARY_PATH="$AMBERHOME/lib:$LD_LIBRARY_PATH"

#Test installation
RUN pmemd -O -h && \
    cpptraj -h

#python
RUN pip3 install numpy pandas matplotlib seaborn mdtraj

WORKDIR /workspace
CMD ["/bin/bash"]
