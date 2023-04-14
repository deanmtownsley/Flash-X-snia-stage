module savelist

export H5_DIR=/nfs/gce/projects/FLASH5/software/hdf5
export VOL_DIR=/nfs/gce/projects/FLASH5/software/vol-async
export ABT_DIR=/nfs/gce/projects/FLASH5/software/argobots

module load autoconf/2 automake/1.16 cmake/3.20 libtool/2.4

pathmungeany /nfs/gce/software/custom/linux-ubuntu20.04-x86_64/anaconda3/rolling/lib 	before LD_LIBRARY_PATH
pathmungeany ${ABT_DIR}/install/lib 	before LD_LIBRARY_PATH
pathmungeany ${VOL_DIR}/src 		before LD_LIBRARY_PATH
pathmungeany ${H5_DIR}/install/lib 	before LD_LIBRARY_PATH

pathmungeany ${H5_DIR}/install/bin 	before PATH

pathmungeany ${ABT_DIR}/install/lib/libabt.so 		before LD_PRELOAD
pathmungeany ${H5_DIR}/install/lib/libhdf5_fortran.so:${VOL_DIR}/src/libh5async.so 		before LD_PRELOAD
pathmungeany ${H5_DIR}/install/lib/libhdf5.so 		before LD_PRELOAD

module list

