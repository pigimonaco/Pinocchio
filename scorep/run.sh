#!/bin/bash

# This script is meant to set the parameters
# to compile and perform code profiling and
# tracing using Score-p tool (if requested)

if [ ${PROFILING} -eq 1 ] && [ ${TRACING} -eq 1 ]
then
    printf "\n\t PROFILING or TRACING exclusively\n"
    exit 1
fi

if [ ${#COMPILER_CC} -eq 0 ]
then
    printf "\n\t Set the COMPILER_CC ...aborting...\n"
    exit 2
fi

if [ ${#COMPILER_CPP} -eq 0 ]
then
    printf "\n\t Set the COMPILER_CPP ...aborting...\n"
    exit 3
fi

if [ ${#SYSTEM} -eq 0 ]
then
    printf "\n\t Set the SYSTEM ...aborting...\n"
    exit 4
fi

if [ ${#EXEC[@]} -eq 0 ]
then
    printf "\n\t EXEC empty ...aborting...\n"
    exit 5
fi

if [ ${#NTASKS[@]} -eq 0 ]
then
    printf "\n\t NTASKS empty ...aborting...\n"
    exit 6
fi

if [ ${#MPI_MAP_BY[@]} -eq 0 ]
then
    printf "\n\t MPI_MAP_BY empty ...aborting...\n"
    exit 7
fi

if [ ${#OMP_THR[@]} -eq 0 ]
then
    printf "\n\t OMP_THR empty ...aborting...\n"
    exit 8
fi

if [ "${OUTPUT_DIRECTORY}" = "${WORKDIR}/scorep" ] || [ "${OUTPUT_DIRECTORY}" = "${WORKDIR}/src" ] || [ "${OUTPUT_DIRECTORY}" = "${WORKDIR}/example" ] || [ "${OUTPUT_DIRECTORY}" = "${WORKDIR}/src/energy" ]
then
    printf "\n\t Invalid output directory: ${OUTPUT_DIRECTORY} ...aborting...\n"
    exit 9
fi

# create the OUTPUT_DIRECTORY
mkdir -p ${OUTPUT_DIRECTORY}

if [ ${PROFILING} -eq 1 ] || [ ${TRACING} -eq 1 ]
then
    DIRNAME_PREFIX=profiling
fi

if [ ${ENERGY} -eq 1 ]
then
    DIRNAME_PREFIX=energy
fi

if [ ${PRODUCTION} -eq 1 ]
then
    DIRNAME_PREFIX=production
fi

# Number of nodes (set automatically by Slurm)
if [ ${#SLURM_NNODES} -eq 0 ]
then
    NODES=1
else
    NODES=${SLURM_NNODES}
fi

# PAPI metrics
PAPI_METRICS=()
# set the output directory
DIR=${OUTPUT_DIRECTORY}/${DIRNAME_PREFIX}_BoxSize_${box}_GridSize_${grid}
OUT_DIR=()

if [ ${PROFILING} -eq 1 ] || [ ${TRACING} -eq 1 ]
then
    if [ ${#PAPI_METRIC} -ne 0 ]
    then
	# get the number of items
	ITEMS=$(echo ${PAPI_METRIC} | gawk -F "," '{print NF}')
	for ((item=0 ; item<${ITEMS} ; item+=2))
	do
	    if [[ ${item} -ne $((ITEMS - 1)) ]]
	    then
		# extract the papi metric
		METRIC=$(echo ${PAPI_METRIC} | gawk -F "," -v start=$((item + 1)) 'BEGIN {stop = start + 1} {print $start","$stop}' | tr -d ' ')
		PAPI_DIR=$(echo ${METRIC} | sed -r 's/,/_/g')
	    else
		METRIC=$(echo ${PAPI_METRIC} | gawk -F "," -v start=$((item + 1)) '{print $start}' | tr -d ' ')
		PAPI_DIR=${METRIC}
	    fi
	    PAPI_METRICS+=( "${METRIC}" )
	    DIR_NAME="${DIR}_${PAPI_DIR}"
	    # make directory
	    mkdir -p ${DIR_NAME}
	    OUT_DIR+=( ${DIR_NAME} )
	done
    else
	# make directory
	mkdir -p ${DIR}
	OUT_DIR+=( ${DIR} )
    fi
fi

if [ ${PRODUCTION} -eq 1 ]
then
    # make directory
    mkdir -p ${DIR}
    OUT_DIR+=( ${DIR} )
fi
