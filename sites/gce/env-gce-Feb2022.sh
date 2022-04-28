# source this file before running on GCE machines with Ubuntu 20.04
module load autoconf/2 automake/1.16 visit/3  cmake/3.20 libtool/2.4 #  git-lfs/2.11.0
export H5_DIR=/nfs/gce/projects/FLASH5/software/hdf5
export LD_LIBRARY_PATH=/nfs/gce/projects/FLASH5/software/hdf5/install/lib:/nfs/gce/software/custom/linux-ubuntu18.04-x86_64/anaconda3/rolling/lib:$LD_LIBRARY_PATH
export PATH=$H5_DIR/install/bin:$PATH
