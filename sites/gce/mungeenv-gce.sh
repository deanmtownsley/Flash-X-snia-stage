module savelist

export H5_DIR=/nfs/gce/projects/FLASH5/software/hdf5

module load autoconf/2.69-tz6eue5 automake/1.16.3-fm5m6qc visit/3.0.0 git-lfs/2.11.0  cmake/3.20.0-vov726r libtool/2.4.6-jdxbjft

pathmungeany /nfs/gce/software/custom/linux-ubuntu18.04-x86_64/anaconda3/rolling/lib 	before LD_LIBRARY_PATH

pathmungeany ${H5_DIR}/install/lib 	before LD_LIBRARY_PATH

pathmungeany ${H5_DIR}/install/bin 	before PATH

module list
