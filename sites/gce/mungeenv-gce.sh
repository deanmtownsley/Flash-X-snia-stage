# Prepare GCE-UB18 machine for compiling and running Flash-X - GCC
module savelist

export H5_DIR=/nfs/gce/projects/FLASH5/software/hdf5

module purge

##pathmungeany /nfs/gce/software/custom/linux-ubuntu18.04-x86_64/anaconda3/rolling/lib 	before LD_LIBRARY_PATH

pathmungeany ${H5_DIR}/install/lib 	before LD_LIBRARY_PATH

pathmungeany ${H5_DIR}/install/bin 	before PATH

module list
