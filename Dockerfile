FROM nvidia/cuda:11.7.1-devel-ubuntu20.04

ARG DEBIAN_FRONTEND=noninteractive

#dependencies
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
    libnetcdf-dev \
    libopenmpi-dev openmpi-bin && \
    rm -rf /var/lib/apt/lists/*

#setup conda
WORKDIR /opt
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh

#set conda environment variables
ENV PATH="/opt/conda/bin:$PATH"
RUN conda init && conda config --set always_yes yes

#activate conda environment (python 3.9)
RUN conda create -n amber-env python=3.9 && \
    echo "conda activate amber-env" >> ~/.bashrc

#get latest ambertools from conda
RUN /opt/conda/bin/conda install -n amber-env -c conda-forge ambertools openmm mdtraj

#set up environment
ENV AMBERHOME="/opt/conda/envs/amber-env"
ENV PATH="$AMBERHOME/bin:$PATH"
ENV LD_LIBRARY_PATH="$AMBERHOME/lib:$LD_LIBRARY_PATH"

#verify installation
RUN /opt/conda/bin/conda run -n amber-env pmemd -O -h && \
    /opt/conda/bin/conda run -n amber-env cpptraj -h

#reset to default workspace
WORKDIR /workspace
CMD ["/bin/bash"]
