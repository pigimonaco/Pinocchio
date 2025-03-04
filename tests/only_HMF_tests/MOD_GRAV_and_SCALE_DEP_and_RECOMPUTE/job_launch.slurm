#!/bin/bash  
#SBATCH --partition=pleiadi
#SBATCH --job-name=Gadget_Ics
#SBATCH --ntasks=36
#SBATCH --ntasks-per-node=36
#SBATCH --cpus-per-task=1  # usually, obviously, in the range[1-10]
#SBATCH --nodes=1
##SBATCH --mem-per-cpu=6gb
#SBATCH -t 0-8:00 # time (D-HH:MM)  
#SBATCH -o slurm.%N.%j.out # STDOUT  
#SBATCH -e slurm.%N.%j.err # STDERR  
#SBATCH --account=pinocchio
 

export MODULE_VERSION=4.0.1 && source /opt/cluster/spack/share/spack/setup-env.sh && module load default-gcc-11.2.0 && module load pfft-1.0.8-alpha-gcc-11.2.0-qy6lqer
mpirun -np 4 ../Pinocchio/src/pinocchio.x parameter_file_standard > log_MOD_GRAV_and_SCALE_DEP_and_RECOMPUTE
