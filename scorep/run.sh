#!/bin/bash

# This script is meant to set the parameters
# to compile and perform code profiling and
# tracing using Score-p tool

###################################### TO EDIT ######################################################
# set the compiler
COMPILER=mpicc

# set the machine name on the src/Makefile accordingly
SYSTEM=LeonardoBoost

# executable basename
BASE=Pinocchio
# executable with/without OMP or Debug support
EXEC=( "${BASE}" "${BASE}OMP" )

# Select the PAPI metrics (empty means do not use PAPI counters)
# PAPI_METRIC=PAPI_DP_OPS,PAPI_L3_TCM
PAPI_METRIC=

# Number of MPI tasks (array)
NTASKS=( 8 )

# MPI mapping list (e.g. socket, numa, core)
MPI_MAP_BY=( numa )

# Number of OMP threads per MPI task (array)
OMP_THR=( 1 )

# Parameter file
PARAMFILE=${WORKDIR}/example/parameter_file_profiling

# Outputs (name of file with required output redshifts)
OUTPUTS=${WORKDIR}/example/outputs_scorep

# Profiling/tracing directory (output directory)
OUT_DIR=${WORKDIR}/profiling_256_par_256_box

####################################### END TO EDIT #################################################
#####################################################################################################
#####################################################################################################

# Number of nodes (set automatically by Slurm)
if [ ${#SLURM_NNODES} -eq 0 ]
then
    NODES=1
else
    NODES=${SLURM_NNODES}
fi

# check
if [ ! -f ${PARAMFILE} ]
then
    printf "\n\t Paramfile: ${PARAMFILE} not found ...aborting...\n"
    exit 0
fi

if [ ! -f ${OUTPUTS} ]
then
    printf "\n\t Outputs: ${OUTPUTS} not found ...aborting...\n"
    exit 1
fi

if [ "${OUT_DIR}" = "${WORKDIR}/scorep" ] || [ "${OUT_DIR}" = "${WORKDIR}/src" ] || [ "${OUT_DIR}" = "${WORKDIR}/example" ]
then
    printf "\n\t Invalid output directory: ${OUT_DIR} ...aborting..."
    exit 2
else
    # create the directory. If already exists then keep files
    mkdir -p ${OUT_DIR}
fi
