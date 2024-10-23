FROM nvcr.io/hpc/gromacs:2023.2

ARG DEBIAN_FRONTEND=noninteractive

## grab an OpenMM installation
RUN apt-get update && apt-get install -y \
    build-essential \
    libssl-dev \
    uuid-dev \
    libgpgme11-dev \
    squashfs-tools \
    libseccomp-dev \
    pkg-config \
    git-all \
    wget \
    libbz2-dev \
    zlib1g-dev \
    python3-dev \
    libffi-dev && \
    mkdir /GW_Python && \
    cd /GW_Python && \
    wget http://www.python.org/ftp/python/3.8.10/Python-3.8.10.tgz && \
    tar -zxvf Python-3.8.10.tgz && \
    cd Python-3.8.10 && \
    ./configure --prefix=/GW_Python && \ 
    cd /GW_Python/Python-3.8.10 && \
    make && \
    make install && \
    /GW_Python/bin/pip3 install OpenMM && \
    #this symlink should hopefully allow gamd-openmm to install to GW_Python
    ln -s /GW_Python/bin/python /usr/bin/env/python && \
    git clone https://github.com/MiaoLab20/gamd-openmm.git /gamd-openmm && \
    ls / && \
    ls ~/ && \
    cd /gamd-openmm && \
    /gamd-openmm/setup.py install
    #/GW_Python/bin/pip3 install torch torchvision torchaudio && \
