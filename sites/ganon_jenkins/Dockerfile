FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update
RUN apt-get install -y python3 python3-pip gfortran build-essential libhdf5-openmpi-dev openmpi-bin pkg-config libopenmpi-dev openmpi-bin libblas-dev liblapack-dev libpnetcdf-dev git python-is-python3
RUN pip3 install numpy matplotlib h5py scipy configparser

# Install amrex here essentially as a cache
RUN git clone https://github.com/AMReX-Codes/amrex.git /amrex
WORKDIR /amrex
RUN git checkout bf267b0954263f9f7d9c340c74bf50efa070f8cc
RUN ./configure --prefix=install/1D --dim=1 && make -j && make install
RUN ./configure --prefix=install/2D --dim=2 && make -j && make install
RUN ./configure --prefix=install/3D --dim=3 && make -j && make install

# The following is needed for the Flash-X-Test script to run. It apparently requires a user string from getpwuid() in python
ENV USER=jenkins
ENV LOGNAME=jenkins
RUN adduser --uid 128 $USER
USER $USER
