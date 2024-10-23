#!/bin/bash

#SBATCH --job-name=PinocchioScorep
#SBATCH --partition=pleiadi
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --output=PinocchioScorep-%j.out
#SBATCH --error=PinocchioScorep-%j.err
#SBATCH --time=08:00:00


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

#------- set the working dir
export WORKDIR=/u/dgoz/beegfs/Pinocchio

#------ switch on/off profiling
PROFILING=1

#------ switch on/off tracing
# Tracing is meant to be run after the profiling
TRACING=0

#------ set the scorep memory if tracing is enable
export SCOREP_TOTAL_MEMORY=11MB # set this value according to estimated memory requirements during profiling analysis

# set the compiler
COMPILER=mpicc

# set the machine name on the src/Makefile accordingly
SYSTEM=pleiadi

# executable basename
BASE=Pinocchio

# executable with/without OMP or Debug support
EXEC=( "${BASE}" "${BASE}OMP" )

# Select the PAPI metrics (empty means do not use PAPI counters)
PAPI_METRIC=PAPI_DP_OPS,PAPI_L3_TCM

# Number of MPI tasks (array)
NTASKS=( 8 )

# MPI mapping list (e.g. socket, numa, core)
MPI_MAP_BY=( numa )

# Number of OMP threads per MPI task (array)
OMP_THR=( 2 4 )

# Parameter file
PARAMFILE=${WORKDIR}/example/parameter_file_profiling

# Outputs (name of file with required output redshifts)
OUTPUTS=${WORKDIR}/example/outputs_scorep

# Profiling/tracing directory (output directory)
OUT_DIR=${WORKDIR}/profiling_256_par_256_box


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

if [ ${PROFILING} -eq 1 ]
then
    printf "\n\t Performing the profiling of Pinocchio \n\n"

    # compile the code with score-p support
    source ${SCOREP_PATH}/compile.sh ${COMPILER}

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
