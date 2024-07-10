#!/bin/bash
#SBATCH --mem-per-cpu 900
#SBATCH -n 2
#SBATCH --cpus-per-task 4
#SBATCH -o stdeo.txt
#SBATCH -p main
#SBATCH --qos main

# above -n and --cpu-per-task are for 4 threads with 2 ranks for a total of 8 threads on 8 cpus

# example queuescript for uahpc
# submit with :
#   sbatch queuescript.sh

export OMP_NUM_THREADS=4
export MESASDK_ROOT=/dmt/townsley/code/mesa/mesasdk
source ${MESASDK_ROOT}/bin/mesasdk_init.sh

export MESA_DIR=/dmt/townsley/code/mesa/mesa_9793


srun ./flash4
