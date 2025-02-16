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

module load openmpi/mlnx/gcc/64/4.1.5rc2

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

export MESASDK_ROOT=/grps2/dmtownsley/common/mesa/mesasdk-23.7.3
source ${MESASDK_ROOT}/bin/mesasdk_init.sh

export MESA_DIR=/grps2/dmtownsley/townsley/code/mesa/mesa-r23.05.1-mesaburn

mpirun ./flashx
