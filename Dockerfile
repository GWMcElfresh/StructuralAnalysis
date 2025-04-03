FROM nvidia/cuda:11.7.1-devel-ubuntu20.04

ARG DEBIAN_FRONTEND=noninteractive

#dependencies
RUN apt-get update && apt-get install -y \
    software-properties-common && \
    add-apt-repository -y ppa:ubuntu-toolchain-r/test && \
    apt-get update && apt-get install -y \
    gcc-11 g++-11 libstdc++6 \
    build-essential \
    libboost-all-dev swig \
    sudo \
    cmake \
    gfortran \
    wget \
    git \
    libnetcdf-dev \
    libopenmpi-dev openmpi-bin \
    libgl1-mesa-glx \
    vim \
    python2.7 \
    tk-dev \
    tcl-dev \
    && rm -rf /var/lib/apt/lists/*

#setup MGLTools (requires python 2.7 - oof)
WORKDIR /opt
RUN wget https://ccsb.scripps.edu/mgltools/download/491/mgltools_x86_64Linux2_1.5.7p1.tar.gz && \
    tar -xvzf mgltools_x86_64Linux2_1.5.7p1.tar.gz && \
    cd mgltools_x86_64Linux2_1.5.7 && \
    yes yes | ./install.sh && \
    rm ../mgltools_x86_64Linux2_1.5.7p1.tar.gz

#pop MGLTools into PATH
ENV PATH="/opt/mgltools_x86_64Linux2_1.5.7/bin:/opt/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24:$PATH"

#setup conda
WORKDIR /opt
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh

#set conda environment variables
ENV PATH="/opt/conda/bin:$PATH"
RUN conda init && conda config --set always_yes yes

#activate conda environment (python 3.9)
RUN conda create -n openmm-env python=3.9.21 && \
    echo "conda activate openmm-env" >> ~/.bashrc

#get latest ambertools from conda
RUN conda install -n openmm-env -c conda-forge openmm ambertools mdtraj numpy scipy pandas matplotlib seaborn mdanalysis seaborn libstdcxx-ng && \
#symlink c++ lib for `GLIBCXX_3.4.30' not found
    sudo ln -sf /usr/lib/x86_64-linux-gnu/libstdc++.so.6 /opt/conda/envs/openmm-env/lib/

#conda init for the vina image
RUN conda init bash

#reload shell
SHELL ["bash", "-c"]

#install a vina env
RUN bash -c "source /opt/conda/etc/profile.d/conda.sh && \
    conda create -n vina python=3.9.21 && \
    conda activate vina && \
    conda config --env --add channels conda-forge && \
    conda install -c conda-forge numpy swig boost-cpp libboost sphinx sphinx_rtd_theme && \
    pip install vina"

RUN bash -c "source /opt/conda/etc/profile.d/conda.sh && \
    conda create -n mordred python=3.9.21 && \
    conda activate mordred && \
    conda config --env --add channels conda-forge && \
    conda install -c rdkit -c mordred-descriptor mordred"

#verify install:
RUN /bin/bash -c 'source activate openmm-env && python3 -c "import openmm; print(openmm.version.version)"'

#reset to default workspace
WORKDIR /workspace
CMD ["/bin/bash"]
