# Dockerfile for setting up Flash-X testing environment.
# Artifacts are cached if this file does not change between commits
#
# This file should be used for one-time setups that do not need 
# to go through changes during every commit. For testing diffs
# in Flash-X source code we use the Jenkinsfile.

# Start from ubuntu base image
FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive

# Upgrade and update distribution
RUN apt-get dist-upgrade
RUN apt-get update

# Installation of compilers, tools, and libraries from ubuntu distribution
RUN apt-get install -y python3 python3-pip gfortran build-essential 
RUN apt-get install -y libhdf5-openmpi-dev openmpi-bin pkg-config libopenmpi-dev openmpi-bin 
RUN apt-get install -y libblas-dev liblapack-dev libpnetcdf-dev git python-is-python3
RUN apt-get install -y curl tar

# Installation of python packages
RUN pip3 install numpy matplotlib h5py scipy configparser
RUN pip3 install git+https://github.com/Flash-X/Flash-X-Test.git@2023.01

# Installation of AMReX from source
RUN git clone https://github.com/Box-Tools/amrex.git /amrex
WORKDIR /amrex
RUN git checkout c68734c
RUN ./configure --prefix=install/1D --dim=1 && make -j && make install
RUN ./configure --prefix=install/2D --dim=2 && make -j && make install
RUN ./configure --prefix=install/3D --dim=3 && make -j && make install

# Install hypre from source
RUN git clone https://github.com/hypre-space/hypre.git /hypre
WORKDIR /hypre
RUN git checkout v2.22.0
RUN cd src && ./configure --enable-shared --enable-fortran --with-MPI --prefix=/hypre/install && make install
ENV LD_LIBRARY_PATH=/hypre/install/lib:$LD_LIBRARY_PATH
ENV DYLD_LIBRARY_PATH=/hypre/install/lib:$DYLD_LIBRARY_PATH

# Create an archive directory to store baselins
RUN mkdir -pv /archive/ganon_jenkins

# Copy benchmarks from cloud archive and unpack them into archive
RUN curl -s -L https://anl.box.com/shared/static/ybxzm1uhgdotxro7q700mahv6gw82ufr | tar xvz -C /archive/ganon_jenkins
RUN curl -s -L https://anl.box.com/shared/static/a4jgr0rxj8d2l3oqfg8zwcgomscsvliq | tar xvz -C /archive/ganon_jenkins
RUN curl -s -L https://anl.box.com/shared/static/gu7qm8xq6pfcpcdp0s8vggumdxwg34fl | tar xvz -C /archive/ganon_jenkins

# Change permissions for the archive directory
RUN chmod -R 777 /archive/ganon_jenkins/* && ls -l /archive/ganon_jenkins

# Run container in non-root mode. This requirement is imposed by Flash-X-Test for security.
ENV USER=jenkins
ENV LOGNAME=jenkins
RUN adduser --uid 128 $USER
USER $USER
