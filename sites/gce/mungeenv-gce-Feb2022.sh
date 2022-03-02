module savelist

export H5_DIR=/nfs/gce/projects/FLASH5/software/hdf5

module load autoconf/2 automake/1.16 cmake/3.20 libtool/2.4
module add anaconda3/rolling


pathmungeany ${H5_DIR}/install/lib 	before LD_LIBRARY_PATH

pathmungeany ${H5_DIR}/install/bin 	before PATH

module list
