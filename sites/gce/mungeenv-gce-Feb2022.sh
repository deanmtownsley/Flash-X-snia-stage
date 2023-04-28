# Prepare GCE-UB20 machine for compiling and running Flash-X - Intel compilers
module savelist

module purge
module add intel/20.4
module add mpich/3.4.2-intel
module add hdf5/1.12.1-mpich-3.4.2-parallel-fortran

module list
