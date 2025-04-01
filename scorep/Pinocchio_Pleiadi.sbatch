#!/bin/bash

#SBATCH --job-name=Pinocchio_Scorep
#SBATCH --partition=pleiadi
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --output=Pinocchio_Scorep-%j.out
#SBATCH --error=Pinocchio_Scorep-%j.err
#SBATCH --time=03:00:00

########################################### MODULES / LIBRARIES #################################################
export MODULE_VERSION=5.0.1
source /opt/cluster/spack/share/spack/setup-env.sh

module purge
module load default-gcc-11.2.0
module load scorep-7.0-gcc-11.2.0-4eap4t5
module load papi-6.0.0.1-gcc-11.2.0-hnhu57s

OTF2=/opt/cluster/spack/opt/spack/linux-centos7-broadwell/gcc-11.2.0/otf2-2.3-wsvfcmjwktwfg3fgkyauk3arsqivj5b2/lib
CUB4W=/opt/cluster/spack/opt/spack/linux-centos7-broadwell/gcc-11.2.0/cubew-4.6-5cfpxt7b3lma5faqiwyep5phk75vqikg/lib/
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${OTF2}:${CUB4W}

############################################## TO EDIT ###########################################################
##################################################################################################################

# export OMP_TARGET_OFFLOAD=mandatory # NOT USED by Pinocchio
export OMP_PROC_BIND=true
export OMP_WAIT_POLICY=ACTIVE

#------- set the working dir
export WORKDIR=/u/dgoz/work/Pinocchio

#-------- swith on/off production
PRODUCTION=0

#-------- swith on/off energy (PMT library is required and code instrumentation)
ENERGY=0  # not working yet

#------ switch on/off profiling
PROFILING=0

#------ switch on/off tracing
# Tracing is meant to be run after the profiling
TRACING=1
#------ set the scorep memory if tracing is enable
export SCOREP_TOTAL_MEMORY=21MB # set this value according to estimated memory requirements during profiling analysis

# set the compilers
COMPILER_CC=mpicc
COMPILER_CPP=mpicxx
# set the machine name on the src/Makefile accordingly
SYSTEM=pleiadi

# executable basename
BASE=Pinocchio

# executable with/without OMP or Debug support
EXEC=( "${BASE}_OMP" "${BASE}" )

# Select the PAPI metrics (empty means do not use PAPI counters)
# string comma-separated
# PAPI_METRIC="PAPI_L1_DCM, PAPI_L2_DCM, PAPI_L3_TCM, PAPI_TLB_DM, PAPI_MEM_SCY, PAPI_TOT_INS, PAPI_TOT_CYC"
# PAPI_METRIC="PAPI_L1_DCM, PAPI_L2_DCM, PAPI_TLB_DM"

# Number of MPI tasks (array)
NTASKS=( 4 8 )

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
OUTPUT_DIRECTORY=${WORKDIR}/example/scorep_no_PAPI

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
