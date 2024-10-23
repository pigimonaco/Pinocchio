#!/bin/bash

if [ $# -lt 1 ]
then
    printf "\n\t Usage: ./compile.sh <compiler (e.g. gcc, nvcc, ecc.)> \n\n"
    exit 1
else
    CC=$1
    if ! command -v ${CC} &>/dev/null
    then
	printf "\n\t ${CC} not found... aborting... \n"
	exit 2
    fi

    if ! command -v scorep &>/dev/null
    then
	printf "\n\t scorep not found... aborting... \n"
	exit 3
    fi
fi

printf "\n\n\t Using the ${CC} compiler \n\n"

cd ${WORKDIR}/src

for EXE in ${EXEC[@]}
do
    DEBUG_FLAG=$(echo ${EXE} | grep "Debug")
    OMP_FLAG=$(echo ${EXE}   | grep "OMP")
    GPU_FLAG=$(echo ${EXE}   | grep "GPU")

    if [ ${#DEBUG_FLAG} -eq 0 ]
    then
	DEBUG_SWITCH=NO
    else
	DEBUG_SWITCH=YES
    fi

    if [ ${#OMP_FLAG} -eq 0 ]
    then
	OMP_SWITCH=NO
    else
	OMP_SWITCH=YES
    fi

    if [ ${#GPU_FLAG} -eq 0 ]
    then
	GPU_SWITCH=NO
    else
	GPU_SWITCH=YES
    fi

    if [ ${PROFILING} -eq 1 ] || [ ${TRACING} -eq 1 ]
    then
	# Scorep compilation
	EXE_SCOREP=${EXE}_Scorep
	SCOREP_SWITCH=YES
    else
	EXE_SCOREP=${EXE}_Prod
	SCOREP_SWITCH=NO
    fi

    printf "\n\t Compiling ${EXE_SCOREP}... \n"
    COMPILE=compile_${EXE_SCOREP}.txt
    make clean && make EXEC=${EXE_SCOREP} COMPILER=${CC} DEBUG=${DEBUG_SWITCH} OMP=${OMP_SWITCH} GPU=${GPU_SWITCH} SYSTYPE=${SYSTEM} SCOREP=${SCOREP_SWITCH} |& tee ${COMPILE}
    file ${EXE_SCOREP} |& tee -a ${COMPILE}
    ldd ${EXE_SCOREP} |& tee -a ${COMPILE}

    for dir in ${OUT_DIR[@]}
    do
	cp ${EXE_SCOREP} ${COMPILE} ${dir}
    done
    make clean
    rm -f ${EXE_SCOREP} ${COMPILE}
done # loop over EXE

cd ${WORKDIR}

printf "\n\t Compilation done \n"
