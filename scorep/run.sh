#!/bin/bash

# This script is meant to set the parameters
# to compile and perform code profiling and
# tracing using Score-p tool

if [ ${PROFILING} -eq 0 ] && [ ${TRACING} -eq 0 ]
then
    printf "\n\t PROFILING or TRACING \n"
    exit 0
fi

if [ ${PROFILING} -eq 1 ] && [ ${TRACING} -eq 1 ]
then
    printf "\n\t PROFILING or TRACING exclusively\n"
    exit 1
fi

if [ ${#COMPILER} -eq 0 ]
then
    printf "\n\t Set the COMPILER ...aborting...\n"
    exit 2
fi

if [ ${#SYSTEM} -eq 0 ]
then
    printf "\n\t Set the SYSTEM ...aborting...\n"
    exit 3
fi

if [ ${#EXEC[@]} -eq 0 ]
then
    printf "\n\t EXEC empty ...aborting...\n"
    exit 4
fi

if [ ${#NTASKS[@]} -eq 0 ]
then
    printf "\n\t NTASKS empty ...aborting...\n"
    exit 5
fi

if [ ${#MPI_MAP_BY[@]} -eq 0 ]
then
    printf "\n\t MPI_MAP_BY empty ...aborting...\n"
    exit 6
fi

if [ ${#OMP_THR[@]} -eq 0 ]
then
    printf "\n\t OMP_THR empty ...aborting...\n"
    exit 7
fi

if [ ! -f ${PARAMFILE} ]
then
    printf "\n\t Paramfile: ${PARAMFILE} not found ...aborting...\n"
    exit 8
fi

if [ ! -f ${OUTPUTS} ]
then
    printf "\n\t Outputs: ${OUTPUTS} not found ...aborting...\n"
    exit 9
fi

if [ "${OUT_DIR}" = "${WORKDIR}/scorep" ] || [ "${OUT_DIR}" = "${WORKDIR}/src" ] || [ "${OUT_DIR}" = "${WORKDIR}/example" ]
then
    printf "\n\t Invalid output directory: ${OUT_DIR} ...aborting..."
    exit 10
else
    # create the directory. If already exists then keep files
    mkdir -p ${OUT_DIR}
fi

# Number of nodes (set automatically by Slurm)
if [ ${#SLURM_NNODES} -eq 0 ]
then
    NODES=1
else
    NODES=${SLURM_NNODES}
fi
