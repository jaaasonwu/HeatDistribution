#!/bin/bash
#SBATCH --job-name=HelloWorld
#SBATCH -p main
#SBATCH --nodes=4
#SBATCH --ntasks=4
#SBATCH --time 0-00:10:00
#SBATCH --cpus-per-task=32

# Run info and srun job launch
srun mpi_heat_distribution 1000 0.001 output.dat
echo 'Done'
