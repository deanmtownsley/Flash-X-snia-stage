FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get dist-upgrade
RUN apt-get update
RUN apt-get install -y python3 python3-pip gfortran build-essential libhdf5-openmpi-dev openmpi-bin pkg-config libopenmpi-dev openmpi-bin libblas-dev liblapack-dev libpnetcdf-dev git python-is-python3
RUN apt-get install -y curl tar
RUN pip3 install numpy matplotlib h5py scipy configparser

# Get archive directories
RUN mkdir -pv /archive/ganon_jenkins

# Copy benchmarks from cloud archive
RUN curl -s -L https://anl.box.com/shared/static/ybxzm1uhgdotxro7q700mahv6gw82ufr | tar xvz -C /archive/ganon_jenkins

# Change permissions
RUN chmod -R 777 /archive/ganon_jenkins/* && ls -l /archive/ganon_jenkins

# Install testing tools
RUN pip3 install git+https://github.com/Flash-X/Flash-X-Test.git@2023.01

# Install amrex here essentially as a cache
RUN git clone https://github.com/Box-Tools/amrex.git /amrex
WORKDIR /amrex
RUN git checkout c68734c
RUN ./configure --prefix=install/1D --dim=1 && make -j && make install
RUN ./configure --prefix=install/2D --dim=2 && make -j && make install
RUN ./configure --prefix=install/3D --dim=3 && make -j && make install

# Install hypre
RUN git clone https://github.com/hypre-space/hypre.git /hypre
WORKDIR /hypre
RUN git checkout v2.22.0
RUN cd src && ./configure --enable-shared --enable-fortran --with-MPI --prefix=/hypre/install && make install

# The following is needed for the Flash-X-Test script to run. It apparently requires a user string from getpwuid() in python
ENV USER=jenkins
ENV LOGNAME=jenkins
RUN adduser --uid 128 $USER
USER $USER
