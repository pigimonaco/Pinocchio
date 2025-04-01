#!/bin/bash

#SBATCH --job-name=Pinocchio-Scorep
#SBATCH --partition=boost_usr_prod
#SBATCH --account=CNHPC_1498509
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:2
#SBATCH --gpu-bind=closest
##SBATCH --exclusive
#SBATCH --output=Pinocchio-Scorep-%j.out
#SBATCH --error=Pinocchio-Scorep-%j.err
#SBATCH --time=00:30:00

module purge
ml gcc/12.2 binutils/2.41 nvhpc/24.3 openmpi/4.1.6--nvhpc--24.3
module load scorep/8.3--openmpi--4.1.6--nvhpc--24.3
# gsl, fftw, pfft have been compiled locally using nvhpc
# src/Makefile has been set up accordingly

# edit the libraries' path
LIB=/leonardo/home/userexternal/dgoz0000/lib/nvhpc-23.11
FTTW_LIB=${LIB}/fftw/fftw-3.3.10/lib
GSL_LIB=${LIB}/gsl/gsl-2.7.1/lib
PFFT_LIB=${LIB}/pfft/pfft/lib
PMT_LIB=/leonardo/home/userexternal/dgoz0000/lib/pmt/local/lib64
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${FTTW_LIB}:${GSL_LIB}:${PFFT_LIB}:${PMT_LIB}

############################################## TO EDIT ###########################################################
##################################################################################################################

# export OMP_TARGET_OFFLOAD=mandatory # NOT USED by Pinocchio
export OMP_PROC_BIND=true
export OMP_WAIT_POLICY=ACTIVE

#------- set the working dir
export WORKDIR=/leonardo/home/userexternal/dgoz0000/FAST_IscrC_SCGPCT/Pinocchio

#-------- swith on/off production
PRODUCTION=1

#-------- swith on/off energy (PMT library is required and code instrumentation)
ENERGY=0  # not working yet

#------ switch on/off profiling
PROFILING=0

#------ switch on/off tracing
# Tracing is meant to be run after the profiling
TRACING=0
#------ set the scorep memory if tracing is enable
export SCOREP_TOTAL_MEMORY=21MB # set this value according to estimated memory requirements during profiling analysis

# set the compilers
COMPILER_CC=nvc
COMPILER_CPP=nvc++
# set the machine name on the src/Makefile accordingly
SYSTEM=LeonardoBoost

# executable basename
BASE=Pinocchio

# executable with/without OMP or Debug support
EXEC=( "${BASE}_OMP" "${BASE}" "${BASE}_OMP_GPU" )

# Select the PAPI metrics (empty means do not use PAPI counters)
# string comma-separated
# PAPI_METRIC= NOT working on LeonardoBoost

# Number of MPI tasks (array)
NTASKS=( 2 4 )

# MPI mapping list (e.g. socket, numa, core)
MPI_MAP_BY=( numa socket )

# Number of OMP threads per MPI task (array)
OMP_THR=( 1 2 )

# box size
BOX=( 128 )

# grid size
GRID=( 128 )

# max memory per node available in Gbytes
MAX_MEM=120

# Parameter file
PARAMFILE_TEMPLATE=${WORKDIR}/example/parameter_file

# Outputs (name of file with required output redshifts)
OUTPUTS=${WORKDIR}/example/outputs

# output directory
OUTPUT_DIRECTORY=${WORKDIR}/example/scorep

####################################################################################################################
####################################################################################################################
#
#
#
#
#
export SCOREP_PATH=${WORKDIR}/scorep

if [ ! -d ${SCOREP_PATH} ]
then
    printf "\n\t ${SCOREP_PATH} does not exist... aborting... \n\n"
    exit 1
fi

# load the run configuration    
source ${SCOREP_PATH}/run.sh

printf "\n\t Running on ${NODES} nodes \n"
printf "\n\t Running Pinocchio on hosts: ${SLURM_NODELIST}"

# get the actual topology
source ${SCOREP_PATH}/topology.sh ${NODES}

if [ ${PRODUCTION} -eq 1 ]
then
    printf "\n\t Performing a production run using Pinocchio \n\n"

    # compile the code without score-p support
    source ${SCOREP_PATH}/compile.sh ${COMPILER_CC} ${COMPILER_CPP}

    # perform the run
    source ${SCOREP_PATH}/production.sh
fi

if [ ${ENERGY} -eq 1 ]
then
    printf "\n\t To do ... \n\n"
fi

if [ ${PROFILING} -eq 1 ]
then
    printf "\n\t Performing the profiling of Pinocchio \n\n"

    # compile the code with score-p support
    source ${SCOREP_PATH}/compile.sh ${COMPILER_CC} ${COMPILER_CPP}

    # perform the profiling
    source ${SCOREP_PATH}/profiling.sh
fi

if [ ${TRACING} -eq 1 ]
then
    printf "\n\t Performing the tracing of Pinocchio \n\n"
    
    source ${SCOREP_PATH}/tracing.sh
fi

printf "\n\t END OF GAME!!! \n\n"

exit 0
