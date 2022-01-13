This directory is for building FLASH on Cray XC40 Theta at ANL. It supports the cray compiler and intel compiler.

For the cray compiler as default, Makefile.h is the same as Makefile.h.cray. Makefile.h.intel is for the intel compiler.

For the cray compiler, before you compile the FLASH, please load the following modules:

module swap PrgEnv-intel PrgEnv-gnu

module load gcc cray-hdf5-parallel cray-parallel-netcdf


For the intel compiler, before you compile the FLASH, please load the following module:

module load gcc cray-hdf5-parallel cray-parallel-netcdf
