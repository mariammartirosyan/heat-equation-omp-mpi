#!/bin/bash
#SBATCH -N 6
#SBATCH --ntasks 96
hostname
srun hostname
mpirun ./a2-mpi --m 1152 --n 1152 --epsilon 0.01 --max-iterations 1000