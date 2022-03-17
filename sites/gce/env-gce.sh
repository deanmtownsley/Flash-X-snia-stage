# source this file before running on GCE machines with Ubuntu 18.04
module load autoconf/2.69-tz6eue5 automake/1.16.3-fm5m6qc visit/3.0.0 git-lfs/2.11.0  cmake/3.20.0-vov726r libtool/2.4.6-jdxbjft
export H5_DIR=/nfs/gce/projects/FLASH5/software/hdf5
export LD_LIBRARY_PATH=/nfs/gce/projects/FLASH5/software/hdf5/install/lib:/nfs/gce/software/custom/linux-ubuntu18.04-x86_64/anaconda3/rolling/lib:$LD_LIBRARY_PATH
export PATH=$H5_DIR/install/bin:$PATH
