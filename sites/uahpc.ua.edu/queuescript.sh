#!/bin/bash
#SBATCH --mem-per-cpu 4096
#SBATCH -n 64
#SBATCH -o stdeo.txt
#SBATCH -p main
#SBATCH --qos main

# example queuescript for uahpc
# submit with :
#   sbatch queuescript.sh

srun ./flash4
