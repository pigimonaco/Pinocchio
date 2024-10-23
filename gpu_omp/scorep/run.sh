#!/bin/bash

# This script is meant to set the parameters
# to compile and perform code profiling and
# tracing using Score-p tool (if requested)

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

if [ ! -f ${PARAMFILE_TEMPLATE} ]
then
    printf "\n\t Paramfile: ${PARAMFILE} not found ...aborting...\n"
    exit 8
fi

if [ ! -f ${OUTPUTS} ]
then
    printf "\n\t Outputs: ${OUTPUTS} not found ...aborting...\n"
    exit 9
fi

if [ "${OUTPUT_DIRECTORY}" = "${WORKDIR}/scorep" ] || [ "${OUTPUT_DIRECTORY}" = "${WORKDIR}/src" ] || [ "${OUTPUT_DIRECTORY}" = "${WORKDIR}/example" ]
then
    printf "\n\t Invalid output directory: ${OUTPUT_DIRECTORY} ...aborting...\n"
    exit 10
fi

# create the OUTPUT_DIRECTORY
mkdir -p ${OUTPUT_DIRECTORY}

if [ ${PROFILING} -eq 1 ] || [ ${TRACING} -eq 1 ]
then
    DIRNAME_PREFIX=profiling
else
    DIRNAME_PREFIX=production
fi

# Number of nodes (set automatically by Slurm)
if [ ${#SLURM_NNODES} -eq 0 ]
then
    NODES=1
else
    NODES=${SLURM_NNODES}
fi

# set the paramfile
PARAMFILE=
# set the output directory
OUT_DIR=
for box in ${BOX[@]}
do
    for grid in ${GRID[@]}
    do 
	DIR=${OUTPUT_DIRECTORY}/${DIRNAME_PREFIX}_BoxSize_${box}_GridSize_${grid} ; mkdir -p ${DIR}
	cp ${OUTPUTS} ${DIR}
	OUT_DIR+=( "${DIR}" )
	NAME=${DIR}/paramfile_BoxSize_${box}_GridSize_${grid}
	gawk '{if ($1=="BoxSize") {$2="'${box}'"} else if ($1=="GridSize") {$2="'${grid}'"} print $0}' ${PARAMFILE_TEMPLATE} > ${NAME}
	PARAMFILE+=( "${NAME}" )
    done # loop over GRID
done # loop over BOX	    

# check
if [ ${#PARAMFILE[@]} -ne ${#OUT_DIR[@]} ]
then
    printf "\n\t Missmatch between paramefile and output directory ...aborting...\n"
    exit 11
fi
