FROM ubuntu:18.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update
RUN apt-get install -y python python3 python3-pip gfortran build-essential libhdf5-openmpi-dev openmpi-bin pkg-config libopenmpi-dev openmpi-bin libblas-dev liblapack-dev libpnetcdf-dev git
RUN pip3 install numpy matplotlib h5py scipy
ENV USER=jenkins
ENV LOGNAME=jenkins
# the following are needed only for test.py, and are not needed in ubuntu 20.04 (though python-is-python3 is needed)
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 10
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
