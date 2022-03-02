export H5_DIR=/nfs/gce/projects/FLASH5/software/hdf5
export VOL_DIR=/nfs/gce/projects/FLASH5/software/vol-async
export ABT_DIR=/nfs/gce/projects/FLASH5/software/argobots

module load autoconf/2 automake/1.16 visit/3  cmake/3.20 libtool/2.4 #  git-lfs/2.11.0
module add anaconda3/rolling

#pathmungeany /nfs/gce/software/custom/linux-ubuntu18.04-x86_64/anaconda3/rolling/lib 	before LD_LIBRARY_PATH
pathmungeany ${H5_DIR}/install/lib 	before LD_LIBRARY_PATH

pathmungeany ${H5_DIR}/install/bin 	before PATH

