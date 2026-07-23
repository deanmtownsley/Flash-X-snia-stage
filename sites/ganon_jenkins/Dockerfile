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
RUN pip3 install git+https://github.com/Flash-X/Flash-X-Test.git@2023.11
RUN chmod 777 /usr/local/lib/python3.8/dist-packages/FlashXTest/backend/FlashTest/ERROR

# Install bittree from source
RUN git clone https://github.com/Flash-X/Bittree.git --branch v2.0 /Bittree
WORKDIR /Bittree
RUN python3 setup.py library --dim 2 --prefix /Bittree/install/2D
RUN cd build && make && make install && cd ..
RUN python3 setup.py library --dim 3 --prefix /Bittree/install/3D
RUN cd build && make && make install && cd ..

# Installation of AMReX from source
RUN git clone https://github.com/AMReX-Codes/amrex.git --branch 23.11 /Amrex
WORKDIR /Amrex
RUN ./configure --prefix=install/1D --dim=1 && make -j && make install
RUN ./configure --prefix=install/2D --dim=2 && make -j && make install
RUN ./configure --prefix=install/3D --dim=3 && make -j && make install

# Install hypre from source
RUN git clone https://github.com/hypre-space/hypre.git --branch v2.22.0 /Hypre
WORKDIR /Hypre
RUN cd src && ./configure --enable-shared --enable-fortran --with-MPI --prefix=/Hypre/install && make install
ENV LD_LIBRARY_PATH=/Hypre/install/lib:$LD_LIBRARY_PATH
ENV DYLD_LIBRARY_PATH=/Hypre/install/lib:$DYLD_LIBRARY_PATH

# Install approximate nearest neighbor (ANN) from source
RUN git clone https://github.com/Box-Tools/ANN-Fortran.git --branch v1.1 /ANN
WORKDIR /ANN
RUN make linux-g++
ENV LD_LIBRARY_PATH=/ANN/lib:$LD_LIBRARY_PATH
ENV DYLD_LIBRARY_PATH=/ANN/lib:$DYLD_LIBRARY_PATH

# Create an archive directory to store baselins
RUN mkdir -pv /archive/ganon_jenkins/0000-00-00
ENV ARCHIVEDIR=/archive/ganon_jenkins/0000-00-00

# Copy benchmarks from cloud archive and unpack them into archive
RUN curl -s -L https://anl.box.com/shared/static/3ic5yepg7ad2m98g6bf4grjynep4aj88 | \
    tar xvz --strip-components=1 -C $ARCHIVEDIR
RUN curl -s -L https://anl.box.com/shared/static/u5lnmqqyqz2d6ms1ramw3i8xggnjbk6d | \
    tar xvz --strip-components=1 -C $ARCHIVEDIR
RUN curl -s -L https://anl.box.com/shared/static/k2dyw9tc1t5myv8dcxoh9op3hdf61el0 | \
    tar xvz --strip-components=1 -C $ARCHIVEDIR
RUN curl -s -L https://anl.box.com/shared/static/o7hvy7hlskmk0db1dwmt7g78v9yqx66x | \
    tar xvz --strip-components=1 -C $ARCHIVEDIR
RUN curl -s -L https://anl.box.com/shared/static/tfjn3r3fq42n1zv5cp6z8b56cjyzdh7j | \
    tar xvz --strip-components=1 -C $ARCHIVEDIR

# For composite tests create a copy of last changed results using the following hack
RUN cp -R $ARCHIVEDIR/Composite* $ARCHIVEDIR/../.

# Change permissions for the archive directory
RUN chmod -R 777 $ARCHIVEDIR && ls -l $ARCHIVEDIR

# Run container in non-root mode. This requirement is imposed by Flash-X-Test for security.
ENV USER=jenkins
ENV LOGNAME=jenkins
RUN adduser --uid 128 $USER
USER $USER
