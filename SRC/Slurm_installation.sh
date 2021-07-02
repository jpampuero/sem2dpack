#!/bin/bash
#SBATCH --partition=x40
#SBATCH --time=1430
#SBATCH -n1
#SBATCH --job-name installation_specfem3D


module purge
module load intel-gnu8-runtime
module load impi
ulimit -s unlimited

make



