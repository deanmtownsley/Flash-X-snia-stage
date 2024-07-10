#!/bin/bash
#SBATCH --mem-per-cpu 4096
#SBATCH -n 64
#SBATCH -o stdeo.txt
#SBATCH -p main
#SBATCH --qos main

# example queuescript for uahpc
# submit with :
#   sbatch queuescript.sh

export MESA_DIR=/your/mesa/dir

export LD_LIBRARY_PATH=$MESA_DIR/utils/hdf5/lib:$LD_LIBRARY_PATH

srun ./flash4
