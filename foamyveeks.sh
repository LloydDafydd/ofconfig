#!/bin/bash

#SBATCH --job-name=foamVTK
#SBATCH --partition=grace
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
# MPI ranks per node (hybrid: ranks * threads = total cores)
#SBATCH --cpus-per-task=1
# OpenMP threads per rank
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH --output=test_case_%j.out
#SBATCH --error=test_case_%j.er

module purge
module load openmpi/5.0.8/5.0.8

cd $HOME/Isambaseball/jobs/results/v085_s2200_m015_a030
foamToVTK
