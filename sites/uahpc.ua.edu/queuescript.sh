#!/bin/bash
#SBATCH --mem-per-cpu 4096
#SBATCH -n 64
#SBATCH -o stdeo.txt
#SBATCH -p main
#SBATCH --qos main

# example queuescript for uahpc
# submit with :
#   sbatch queuescript.sh

module load openmpi/mlnx/gcc/64/4.1.5rc2
module load hdf5/1.12.1

mpirun ./flashx
