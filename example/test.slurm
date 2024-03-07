#!/bin/bash
#SBATCH --nodes=32
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH -p boost_usr_prod    
#SBATCH --account=IscrC_SCGPCT
#SBATCH -J GPU_Test_Pinocchio
#SBATCH --time=10:00:00
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --output %x-slurm_job-%j.out
#SBATCH --error  %x-slurm_job-%j.err

module load gcc/11.3.0 && module load gsl/2.7.1--gcc--11.3.0 && module load openmpi/4.1.4--gcc--11.3.0-cuda-11.8 && module load nvhpc

export OMPI_MCA_btl=^uct

export OMP_PLACES=sockets
export OMP_PROC_BIND=master

mpirun -np 4 ../src/pinocchio_GPU.x parameter_file > log