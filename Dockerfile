FROM nvidia/cuda:11.7.1-devel-ubuntu20.04

ARG DEBIAN_FRONTEND=noninteractive

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

#add conda
# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /miniconda.sh && \
    chmod +x /miniconda.sh && \
    /miniconda.sh -b -p /opt/conda && \
    rm /miniconda.sh

#add conda to PATH
ENV PATH=/opt/conda/bin:$PATH

#AmberTools24
WORKDIR /opt
RUN conda install -c conda-forge ambertools

WORKDIR /opt/amber23
RUN ./configure gnu && \
    make install -j$(nproc)

ENV AMBERHOME=/opt/amber23
ENV PATH="$AMBERHOME/bin:$PATH"
ENV LD_LIBRARY_PATH="$AMBERHOME/lib:$LD_LIBRARY_PATH"

#Test installation
RUN pmemd -O -h && \
    cpptraj -h

#python
RUN pip3 install numpy pandas matplotlib seaborn mdtraj

WORKDIR /workspace
CMD ["/bin/bash"]
